"""models.py -- Modelos de circuito equivalente (ECM) para células Li-ion.

Etapa E3 do plano de trabalho (PSI3591): Rint e Thevenin 1RC, identificação
de parâmetros e métricas de tensão.

Convenções (devem ser as mesmas em todo o projeto e na especificação D3)
-----------------------------------------------------------------------
* Corrente I > 0 na DESCARGA e I < 0 na CARGA (Plett 2004; Xing et al. 2014).
* SoC em fração [0, 1]; capacidade em Ah; tempo em s; tensão em V; R em ohm.
* Amostragem uniforme com período dt. A discretização assume corrente
  constante entre amostras (Plett 2004, Parte 2).
* O índice k do vetor de saída corresponde ao instante k*dt. Os estados
  evoluem como  x[k+1] = f(x[k], i[k])  e a tensão em k usa i[k] e x[k].

Equações implementadas
----------------------
SoC (Plett 2004, Parte 2, eq. 4)::

    z[k+1] = z[k] - (eta_i * dt / (3600 * Q)) * i[k]
    eta_i = 1 na descarga, eta <= 1 na carga

Rint (Hu et al. 2012, modelo 2; Xing et al. 2014, eq. 1)::

    V[k] = OCV(z[k]) - R0 * i[k]

Thevenin 1RC (Hu et al. 2012, modelo 7), com tau1 = R1*C1::

    U1[k+1] = exp(-dt/tau1) * U1[k] + R1 * (1 - exp(-dt/tau1)) * i[k]
    V[k]    = OCV(z[k]) - R0 * i[k] - U1[k]

Referências
-----------
HU, LI, PENG. J. Power Sources 198 (2012) 359-367.
PLETT. J. Power Sources 134 (2004) 252-261, 262-276, 277-292.
XING et al. Applied Energy 113 (2014) 106-115.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import asdict, dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter

__all__ = [
    "OCVTable", "RintParams", "RintModel", "TheveninParams", "Thevenin1RC",
    "SimResult", "FitResult", "sampling_period", "integrate_soc",
    "rc_discretization", "fit_rint", "fit_thevenin_1rc",
    "mae", "rmse", "max_abs_error", "error_summary", "evaluate_voltage",
    "save_params", "load_params",
]


# --------------------------------------------------------------------------
# Utilidades de validação
# --------------------------------------------------------------------------
def _as_1d(x, name: str) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"'{name}' deve ser um vetor 1-D não vazio.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"'{name}' contém NaN ou Inf.")
    return arr


def _same_length(**arrays: np.ndarray) -> None:
    lengths = {k: len(v) for k, v in arrays.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Vetores com tamanhos diferentes: {lengths}")


def sampling_period(t, rtol: float = 1e-3) -> float:
    """Período de amostragem dt [s] de um vetor de tempo uniforme.

    Levanta ValueError se o tempo não for crescente ou não for uniforme
    (reamostre os dados antes: os modelos assumem dt constante).
    """
    t = _as_1d(t, "t")
    if t.size < 2:
        raise ValueError("São necessárias ao menos 2 amostras.")
    d = np.diff(t)
    if np.any(d <= 0):
        raise ValueError("O vetor de tempo deve ser estritamente crescente.")
    dt = float(np.median(d))
    if np.max(np.abs(d - dt)) > rtol * dt:
        raise ValueError(
            "Amostragem não uniforme (desvio > "
            f"{rtol:.1e} do dt mediano = {dt:g} s). Reamostre os dados."
        )
    return dt


# --------------------------------------------------------------------------
# Curva OCV-SoC
# --------------------------------------------------------------------------
class OCVTable:
    """Curva OCV(SoC) em tabela com interpolação linear.

    Plett (2004, Parte 2) usa uma tabela de 200 pontos com interpolação
    linear. Fora do intervalo da tabela o valor é mantido constante (o SoC
    simulado pode sair de [0, 1] se a capacidade/estado inicial estiverem
    errados; verifique isso nos resultados).

    Parameters
    ----------
    soc, ocv : vetores de mesmo tamanho. SoC em fração.
    monotonic : 'raise' (padrão) levanta erro se OCV não for não-decrescente;
        'fix' aplica um máximo acumulado (útil para dados ruidosos, mas
        registre que a correção foi feita).
    """

    def __init__(self, soc, ocv, monotonic: str = "raise"):
        soc = _as_1d(soc, "soc")
        ocv = _as_1d(ocv, "ocv")
        _same_length(soc=soc, ocv=ocv)
        order = np.argsort(soc)
        soc, ocv = soc[order], ocv[order]
        if np.any(np.diff(soc) <= 0):
            raise ValueError("Valores de SoC devem ser distintos.")
        if np.any(np.diff(ocv) < 0):
            if monotonic == "fix":
                ocv = np.maximum.accumulate(ocv)
            else:
                raise ValueError(
                    "OCV não é monotônica em SoC. Suavize os dados ou use "
                    "monotonic='fix' (e documente a correção)."
                )
        self.soc = soc
        self.ocv = ocv

    def __call__(self, soc) -> np.ndarray:
        return np.interp(np.asarray(soc, dtype=float), self.soc, self.ocv)

    def slope(self, soc) -> np.ndarray:
        """dOCV/dSoC [V por unidade de SoC] (derivada da curva linear por
        partes). Necessária para a Jacobiana de observação do EKF."""
        z = np.asarray(soc, dtype=float)
        idx = np.clip(np.searchsorted(self.soc, z, side="right") - 1,
                      0, len(self.soc) - 2)
        return (self.ocv[idx + 1] - self.ocv[idx]) / (
            self.soc[idx + 1] - self.soc[idx])

    @classmethod
    def from_charge_discharge(cls, soc_chg, v_chg, soc_dis, v_dis,
                              n_points: int = 200) -> "OCVTable":
        """Curva OCV como média das curvas de carga e descarga em baixa
        corrente (Plett 2004, Parte 2; Xing et al. 2014), reamostrada em
        ``n_points`` pontos na faixa de SoC comum às duas curvas. A histerese
        é ignorada, como nos dois trabalhos."""
        sc, vc = _as_1d(soc_chg, "soc_chg"), _as_1d(v_chg, "v_chg")
        sd, vd = _as_1d(soc_dis, "soc_dis"), _as_1d(v_dis, "v_dis")
        _same_length(soc_chg=sc, v_chg=vc)
        _same_length(soc_dis=sd, v_dis=vd)
        oc, od = np.argsort(sc), np.argsort(sd)
        sc, vc, sd, vd = sc[oc], vc[oc], sd[od], vd[od]
        lo, hi = max(sc.min(), sd.min()), min(sc.max(), sd.max())
        if hi <= lo:
            raise ValueError("Curvas de carga e descarga sem faixa de SoC comum.")
        grid = np.linspace(lo, hi, n_points)
        ocv = 0.5 * (np.interp(grid, sc, vc) + np.interp(grid, sd, vd))
        return cls(grid, ocv, monotonic="fix")

    def to_dict(self) -> dict:
        return {"soc": self.soc.tolist(), "ocv": self.ocv.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> "OCVTable":
        return cls(d["soc"], d["ocv"])


# --------------------------------------------------------------------------
# SoC por integração de corrente
# --------------------------------------------------------------------------
def integrate_soc(i, dt: float, soc0: float, capacity_ah: float,
                  eta: float = 1.0) -> np.ndarray:
    """SoC por contagem de Coulomb: z[k+1] = z[k] - eta_i*dt*i[k]/(3600*Q).

    ``eta`` (eficiência coulombiana) só é aplicada na carga (i < 0).
    Retorna vetor com o mesmo tamanho de ``i``; ``soc[0] = soc0``.
    """
    i = _as_1d(i, "i")
    if capacity_ah <= 0 or dt <= 0:
        raise ValueError("capacity_ah e dt devem ser positivos.")
    eta_i = np.where(i < 0, eta, 1.0)
    dz = eta_i * i * dt / (3600.0 * capacity_ah)
    soc = np.empty_like(i)
    soc[0] = soc0
    soc[1:] = soc0 - np.cumsum(dz[:-1])
    return soc


# --------------------------------------------------------------------------
# Resultados de simulação
# --------------------------------------------------------------------------
@dataclass
class SimResult:
    """Saída de ``simulate``: tensão, SoC e (se houver) tensão do ramo RC."""
    v: np.ndarray
    soc: np.ndarray
    v1: Optional[np.ndarray] = None


# --------------------------------------------------------------------------
# Modelo Rint
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class RintParams:
    """r0 na descarga; r0_charge opcional (Hu et al. 2012 distinguem R0+ e
    R0-). Se None, usa-se r0 também na carga."""
    r0: float
    r0_charge: Optional[float] = None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "RintParams":
        return cls(**d)


class RintModel:
    """V[k] = OCV(z[k]) - R0 * i[k]  (sem estado dinâmico)."""
    name = "rint"

    def __init__(self, params: RintParams, ocv: OCVTable):
        self.params = params
        self.ocv = ocv

    def _r0(self, i: np.ndarray) -> np.ndarray:
        r_chg = self.params.r0 if self.params.r0_charge is None \
            else self.params.r0_charge
        return np.where(i < 0, r_chg, self.params.r0)

    def voltage(self, i, soc, dt: Optional[float] = None,
                v1_0: float = 0.0) -> np.ndarray:
        """Tensão prevista dado um SoC (por exemplo, o SoC de referência).
        ``dt`` e ``v1_0`` são ignorados (existem para manter a mesma
        interface do 1RC)."""
        i = _as_1d(i, "i")
        soc = _as_1d(soc, "soc")
        _same_length(i=i, soc=soc)
        return self.ocv(soc) - self._r0(i) * i

    def simulate(self, i, dt: float, soc0: float, capacity_ah: float,
                 eta: float = 1.0) -> SimResult:
        soc = integrate_soc(i, dt, soc0, capacity_ah, eta)
        return SimResult(v=self.voltage(i, soc), soc=soc)


# --------------------------------------------------------------------------
# Modelo Thevenin 1RC
# --------------------------------------------------------------------------
def rc_discretization(r1: float, tau1: float, dt: float,
                      method: str = "exact") -> tuple[float, float]:
    """Coeficientes (a, b) de U1[k+1] = a*U1[k] + b*i[k].

    'exact' (zero-order hold; Hu et al. 2012): a = exp(-dt/tau1),
    b = R1*(1 - a). Vale para qualquer dt se i for constante entre amostras.
    'euler' (Euler direto): a = 1 - dt/tau1, b = R1*dt/tau1. Só é adequado
    se dt << tau1; serve para estudar o erro de discretização.
    """
    if tau1 <= 0 or dt <= 0:
        raise ValueError("tau1 e dt devem ser positivos.")
    if method == "exact":
        a = float(np.exp(-dt / tau1))
        return a, r1 * (1.0 - a)
    if method == "euler":
        return 1.0 - dt / tau1, r1 * dt / tau1
    raise ValueError("method deve ser 'exact' ou 'euler'.")


def _rc_voltage(i: np.ndarray, a: float, b: float, v1_0: float) -> np.ndarray:
    """U1 por filtro IIR: U1[0] = v1_0 e U1[k+1] = a*U1[k] + b*i[k].

    Com b0 = 0, a saída em k=0 é o estado inicial zi[0] (transposed DF-II
    do ``scipy.signal.lfilter``), então zi = [v1_0]."""
    y, _ = lfilter([0.0, b], [1.0, -a], i, zi=[v1_0])
    return y


@dataclass(frozen=True)
class TheveninParams:
    """r0 [ohm], r1 [ohm], tau1 = R1*C1 [s]."""
    r0: float
    r1: float
    tau1: float

    @property
    def c1(self) -> float:
        """C1 = tau1 / R1 [F]."""
        return self.tau1 / self.r1

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "TheveninParams":
        return cls(**d)


class Thevenin1RC:
    """Thevenin 1RC com parâmetros constantes (versão inicial do plano:
    parâmetros constantes ou tabelados, desde que testado e documentado)."""
    name = "thevenin_1rc"

    def __init__(self, params: TheveninParams, ocv: OCVTable,
                 discretization: str = "exact"):
        self.params = params
        self.ocv = ocv
        self.discretization = discretization

    def polarization(self, i, dt: float, v1_0: float = 0.0) -> np.ndarray:
        """Tensão U1[k] do ramo RC para uma sequência de corrente."""
        i = _as_1d(i, "i")
        a, b = rc_discretization(self.params.r1, self.params.tau1, dt,
                                 self.discretization)
        return _rc_voltage(i, a, b, v1_0)

    def voltage(self, i, soc, dt: float, v1_0: float = 0.0) -> np.ndarray:
        """V[k] = OCV(z[k]) - R0*i[k] - U1[k], dado um SoC."""
        i = _as_1d(i, "i")
        soc = _as_1d(soc, "soc")
        _same_length(i=i, soc=soc)
        u1 = self.polarization(i, dt, v1_0)
        return self.ocv(soc) - self.params.r0 * i - u1

    def simulate(self, i, dt: float, soc0: float, capacity_ah: float,
                 eta: float = 1.0, v1_0: float = 0.0) -> SimResult:
        soc = integrate_soc(i, dt, soc0, capacity_ah, eta)
        u1 = self.polarization(i, dt, v1_0)
        v = self.ocv(soc) - self.params.r0 * np.asarray(i, float) - u1
        return SimResult(v=v, soc=soc, v1=u1)


# --------------------------------------------------------------------------
# Métricas (Xing et al. 2014, eqs. 3 e 4; plano, seção 20)
# --------------------------------------------------------------------------
def _err(reference, estimate) -> np.ndarray:
    r, e = _as_1d(reference, "reference"), _as_1d(estimate, "estimate")
    _same_length(reference=r, estimate=e)
    return r - e


def mae(reference, estimate) -> float:
    """Erro médio absoluto."""
    return float(np.mean(np.abs(_err(reference, estimate))))


def rmse(reference, estimate) -> float:
    """Raiz do erro quadrático médio."""
    return float(np.sqrt(np.mean(_err(reference, estimate) ** 2)))


def max_abs_error(reference, estimate) -> float:
    """Maior desvio absoluto."""
    return float(np.max(np.abs(_err(reference, estimate))))


def error_summary(reference, estimate) -> dict:
    return {"mae": mae(reference, estimate),
            "rmse": rmse(reference, estimate),
            "max_abs": max_abs_error(reference, estimate)}


def evaluate_voltage(model, i, v, soc, dt: float, v1_0: float = 0.0,
                     skip: int = 0) -> dict:
    """Avalia um modelo (Rint ou 1RC) contra a tensão medida.

    Usa o SoC de referência ``soc`` (não o simulado), isolando o erro do
    modelo de tensão do erro de contagem de Coulomb. ``skip`` descarta as
    primeiras amostras (transitório da condição inicial de U1).
    Retorna MAE/RMSE/erro máximo em V e o resíduo ``v - v_prev`` no tempo
    (para análise de estrutura temporal, como pede o plano).
    """
    v = _as_1d(v, "v")
    v_pred = model.voltage(i, soc, dt, v1_0)
    res = v - v_pred
    s = slice(skip, None)
    out = error_summary(v[s], v_pred[s])
    out["residual"] = res
    return out


# --------------------------------------------------------------------------
# Identificação de parâmetros
# --------------------------------------------------------------------------
@dataclass
class FitResult:
    """Resultado de uma identificação."""
    params: object
    rmse_v: float
    mae_v: float
    n_samples: int
    success: bool = True
    message: str = ""
    at_bounds: list = field(default_factory=list)
    param_std: Optional[dict] = None
    start_costs: Optional[list] = None


def fit_rint(i, v, soc, ocv: OCVTable, directional: bool = False) -> FitResult:
    """Identifica R0 do Rint por mínimos quadrados lineares em forma fechada
    (Plett 2004, Parte 2; Xing et al. 2014).

    Como V - OCV(z) = -R0*i, resolve-se ``lstsq`` em y = v - OCV(soc).
    Com ``directional=True`` estimam-se R0 (descarga, i>0) e R0 (carga, i<0)
    separadamente, como R0+ e R0- em Hu et al. (2012).
    """
    i, v, soc = _as_1d(i, "i"), _as_1d(v, "v"), _as_1d(soc, "soc")
    _same_length(i=i, v=v, soc=soc)
    y = v - ocv(soc)
    if directional:
        cols = [np.where(i > 0, i, 0.0), np.where(i < 0, i, 0.0)]
        if not (np.any(cols[0] != 0) and np.any(cols[1] != 0)):
            raise ValueError("directional=True exige amostras de carga e "
                             "de descarga nos dados.")
        H = -np.column_stack(cols)
    else:
        H = -i[:, None]
    theta, *_ = np.linalg.lstsq(H, y, rcond=None)
    if np.any(theta <= 0):
        warnings.warn("R0 identificado não positivo: verifique a convenção "
                      "de sinal da corrente, a curva OCV e o SoC.")
    params = RintParams(float(theta[0]),
                        float(theta[1]) if directional else None)
    model = RintModel(params, ocv)
    pred = model.voltage(i, soc)
    return FitResult(params=params, rmse_v=rmse(v, pred), mae_v=mae(v, pred),
                     n_samples=len(i))


def fit_thevenin_1rc(i, v, soc, ocv: OCVTable, dt: float,
                     bounds: Optional[tuple] = None,
                     tau_starts=None, v1_0: float = 0.0,
                     discretization: str = "exact") -> FitResult:
    """Identifica (R0, R1, tau1) do Thevenin 1RC minimizando o erro de tensão
    com ``scipy.optimize.least_squares`` (método 'trf', com limites).

    Hu et al. (2012) usaram otimização global por haver mínimos locais; aqui
    usa-se *multi-start* em tau1 (``tau_starts``) e fica-se com a melhor
    solução. Os limites padrão são R0, R1 em [1e-4, 0.5] ohm e tau1 em
    [dt, 2000] s: AJUSTE-OS à célula e registre-os (parâmetros colados em um
    limite sinalizam problema de identificabilidade ou de modelo).

    Assume U1[0] = ``v1_0`` (por padrão 0: ensaio começando em repouso) e
    usa o SoC fornecido (referência), não o simulado.

    ``param_std`` é uma estimativa aproximada (Jacobiana no ótimo) que supõe
    resíduos independentes; como resíduos de ECM são autocorrelacionados, ela
    é otimista -- use apenas como indicador relativo de identificabilidade.
    """
    i, v, soc = _as_1d(i, "i"), _as_1d(v, "v"), _as_1d(soc, "soc")
    _same_length(i=i, v=v, soc=soc)
    if dt <= 0:
        raise ValueError("dt deve ser positivo.")

    if bounds is None:
        lo = np.array([1e-4, 1e-4, dt])
        hi = np.array([0.5, 0.5, 2000.0])
    else:
        lo, hi = (np.asarray(b, float) for b in bounds)
        if lo.shape != (3,) or hi.shape != (3,) or np.any(lo >= hi):
            raise ValueError("bounds deve ser (lo, hi), cada um com 3 valores "
                             "[r0, r1, tau1] e lo < hi.")
    if tau_starts is None:
        tau_starts = np.geomspace(max(10 * dt, lo[2]), min(1000.0, hi[2]), 6)
    ov = ocv(soc)

    def predict(theta):
        r0, r1, tau1 = theta
        a, b = rc_discretization(r1, tau1, dt, discretization)
        return ov - r0 * i - _rc_voltage(i, a, b, v1_0)

    def resid(theta):
        return v - predict(theta)

    # chute inicial de R0 pelo Rint (mesma convenção de sinal)
    r0_init = fit_rint(i, v, soc, ocv).params.r0 if np.any(i != 0) else lo[0]
    margin = 1e-6 * (hi - lo)

    best, costs = None, []
    for tau0 in np.atleast_1d(tau_starts):
        x0 = np.clip([max(0.5 * r0_init, lo[0]), max(0.5 * r0_init, lo[1]),
                      tau0], lo + margin, hi - margin)
        try:
            sol = least_squares(resid, x0, bounds=(lo, hi), method="trf",
                                x_scale="jac")
        except ValueError as exc:  # falha numérica em um ponto inicial
            costs.append(float("nan"))
            warnings.warn(f"Início tau0={tau0:g} s falhou: {exc}")
            continue
        costs.append(float(sol.cost))
        if best is None or sol.cost < best.cost:
            best = sol
    if best is None:
        raise RuntimeError("Todas as inicializações falharam.")

    r0, r1, tau1 = best.x
    params = TheveninParams(float(r0), float(r1), float(tau1))
    pred = predict(best.x)
    names = ["r0", "r1", "tau1"]
    tol = 1e-3 * (hi - lo)
    at_bounds = [n for n, x, l, h in zip(names, best.x, lo, hi)
                 if x - l < tol[names.index(n)] or h - x < tol[names.index(n)]]
    if at_bounds:
        warnings.warn(f"Parâmetros nos limites: {at_bounds}. Revise limites, "
                      "dados de excitação ou estrutura do modelo.")
    std = None
    dof = len(i) - 3
    if dof > 0:
        try:
            s2 = 2.0 * best.cost / dof
            cov = s2 * np.linalg.pinv(best.jac.T @ best.jac)
            std = dict(zip(names, np.sqrt(np.diag(cov)).tolist()))
        except np.linalg.LinAlgError:
            std = None
    return FitResult(params=params, rmse_v=rmse(v, pred), mae_v=mae(v, pred),
                     n_samples=len(i), success=bool(best.success),
                     message=str(best.message), at_bounds=at_bounds,
                     param_std=std, start_costs=costs)


# --------------------------------------------------------------------------
# Persistência de parâmetros (reprodutibilidade)
# --------------------------------------------------------------------------
def save_params(path: str, model_name: str, params, meta: Optional[dict] = None
                ) -> None:
    """Salva parâmetros identificados em JSON, com metadados do ensaio
    (perfil, temperatura, célula, limites, versão do código...)."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"model": model_name, "params": params.to_dict(),
                   "meta": meta or {}}, f, indent=2, ensure_ascii=False)


def load_params(path: str):
    """Lê um JSON de ``save_params``. Retorna (params, meta)."""
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    cls = {"rint": RintParams, "thevenin_1rc": TheveninParams}[d["model"]]
    return cls.from_dict(d["params"]), d.get("meta", {})
