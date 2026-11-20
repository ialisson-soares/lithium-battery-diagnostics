"""Testes de models.py (pytest). Executar com:  pytest -q

Os testes usam apenas dados sintéticos e soluções analíticas, de modo que
qualquer falha indica erro de implementação (e não de dados).
"""
import numpy as np
import pytest
from scipy.signal import cont2discrete

from models import (OCVTable, RintParams, RintModel, TheveninParams,
                    Thevenin1RC, integrate_soc, rc_discretization,
                    fit_rint, fit_thevenin_1rc, mae, rmse, max_abs_error,
                    sampling_period, evaluate_voltage, save_params,
                    load_params)


def make_ocv() -> OCVTable:
    """OCV sintética monotônica, com platô (tipo LiFePO4)."""
    z = np.linspace(0, 1, 201)
    v = (3.2 + 0.15 * z - 0.2 * np.exp(-z / 0.05)
         + 0.25 * np.exp(-(1 - z) / 0.03))
    return OCVTable(z, v)


def make_current(n=6000, dt=1.0, seed=0):
    """Corrente de excitação: pulsos de descarga/carga com repousos."""
    rng = np.random.default_rng(seed)
    i = np.zeros(n)
    k = 0
    while k < n:
        dur = int(rng.integers(30, 200))
        amp = rng.choice([0.0, 0.5, 1.0, 2.0, -0.5, -1.0, 3.0])
        i[k:k + dur] = amp
        k += dur
    return i


# ---------------------------------------------------------------- utilidades
def test_sampling_period_uniform_and_nonuniform():
    assert sampling_period(np.arange(0, 10, 0.5)) == pytest.approx(0.5)
    t = np.array([0, 1, 2, 3.5, 4.5])
    with pytest.raises(ValueError):
        sampling_period(t)
    with pytest.raises(ValueError):
        sampling_period([0, 1, 1, 2])


# --------------------------------------------------------------------- OCV
def test_ocv_interpolation_and_slope():
    t = OCVTable([0, 0.5, 1.0], [3.0, 3.2, 4.0])
    assert t(0.25) == pytest.approx(3.1)
    assert t(-1) == pytest.approx(3.0)      # fora da faixa: constante
    assert t.slope(0.25) == pytest.approx(0.4)
    assert t.slope(0.75) == pytest.approx(1.6)


def test_ocv_rejects_non_monotonic_unless_fixed():
    with pytest.raises(ValueError):
        OCVTable([0, 0.5, 1], [3.0, 3.3, 3.2])
    t = OCVTable([0, 0.5, 1], [3.0, 3.3, 3.2], monotonic="fix")
    assert np.all(np.diff(t.ocv) >= 0)


def test_ocv_from_charge_discharge_is_average():
    z = np.linspace(0, 1, 50)
    t = OCVTable.from_charge_discharge(z, 3.0 + z + 0.01, z, 3.0 + z - 0.01)
    assert np.allclose(t(z), 3.0 + z, atol=1e-9)


# --------------------------------------------------------------- Coulomb
def test_integrate_soc_one_hour_at_1c():
    dt, q = 1.0, 2.0
    i = np.full(3600, q)                    # 1C de descarga
    soc = integrate_soc(i, dt, 1.0, q)
    assert soc[0] == 1.0
    # soc[k] usa i[0..k-1]  ->  após 3599 passos restam 1/3600
    assert soc[-1] == pytest.approx(1.0 / 3600.0, abs=1e-9)


def test_integrate_soc_charge_efficiency():
    soc = integrate_soc([-1.0, -1.0, -1.0], 3600.0, 0.5, 1.0, eta=0.9)
    assert soc[1] == pytest.approx(0.5 + 0.9)  # eta só na carga


# ------------------------------------------------------------ Discretização
def test_rc_discretization_matches_cont2discrete_zoh():
    r1, tau, dt = 0.03, 40.0, 1.0
    a_exact, b_exact = rc_discretization(r1, tau, dt, "exact")
    # dU1/dt = -U1/tau + (R1/tau) i
    Ad, Bd, *_ = cont2discrete((np.array([[-1 / tau]]), np.array([[r1 / tau]]),
                                np.array([[1.0]]), np.array([[0.0]])),
                               dt, method="zoh")
    assert a_exact == pytest.approx(Ad[0, 0], rel=1e-12)
    assert b_exact == pytest.approx(Bd[0, 0], rel=1e-12)


def test_euler_approaches_exact_for_small_dt():
    ae, be = rc_discretization(0.03, 100.0, 0.01, "exact")
    au, bu = rc_discretization(0.03, 100.0, 0.01, "euler")
    assert au == pytest.approx(ae, abs=1e-7)
    assert bu == pytest.approx(be, rel=1e-3)


# ----------------------------------------------------------------- Rint
def test_rint_voltage_and_direction_dependent_r0():
    ocv = make_ocv()
    m = RintModel(RintParams(0.05, 0.08), ocv)
    v = m.voltage([1.0, -1.0], [0.5, 0.5])
    assert v[0] == pytest.approx(ocv(0.5) - 0.05)
    assert v[1] == pytest.approx(ocv(0.5) + 0.08)


def test_fit_rint_recovers_r0_noiseless_and_noisy():
    ocv, dt = make_ocv(), 1.0
    i = make_current()
    soc = integrate_soc(i, dt, 0.95, 5.0)
    truth = RintModel(RintParams(0.045), ocv)
    v = truth.voltage(i, soc)
    assert fit_rint(i, v, soc, ocv).params.r0 == pytest.approx(0.045, rel=1e-9)
    noisy = v + np.random.default_rng(1).normal(0, 1e-3, v.size)
    assert fit_rint(i, noisy, soc, ocv).params.r0 == pytest.approx(0.045, rel=0.02)


def test_fit_rint_directional():
    ocv = make_ocv()
    i = make_current()
    soc = integrate_soc(i, 1.0, 0.95, 5.0)
    v = RintModel(RintParams(0.04, 0.07), ocv).voltage(i, soc)
    p = fit_rint(i, v, soc, ocv, directional=True).params
    assert p.r0 == pytest.approx(0.04, rel=1e-9)
    assert p.r0_charge == pytest.approx(0.07, rel=1e-9)


# --------------------------------------------------------------- 1RC
def test_thevenin_step_response_analytic():
    r1, tau, dt, i0 = 0.03, 40.0, 1.0, 2.0
    ocv = make_ocv()
    m = Thevenin1RC(TheveninParams(0.04, r1, tau), ocv)
    i = np.full(300, i0)
    u1 = m.polarization(i, dt)
    k = np.arange(300)
    expected = r1 * i0 * (1 - np.exp(-k * dt / tau))   # U1[0] = 0
    assert np.allclose(u1, expected, atol=1e-12)


def test_thevenin_initial_state_decays_at_rest():
    m = Thevenin1RC(TheveninParams(0.04, 0.03, 40.0), make_ocv())
    u1 = m.polarization(np.zeros(100), 1.0, v1_0=0.05)
    assert u1[0] == pytest.approx(0.05)
    assert u1[50] == pytest.approx(0.05 * np.exp(-50 / 40.0), rel=1e-12)


def test_thevenin_matches_naive_loop():
    """lfilter deve coincidir com a recursão explícita do artigo (Hu)."""
    p = TheveninParams(0.04, 0.03, 35.0)
    ocv = make_ocv()
    m = Thevenin1RC(p, ocv)
    i = make_current(500)
    sim = m.simulate(i, 1.0, 0.9, 5.0)
    u = np.zeros(500)
    a, b = rc_discretization(p.r1, p.tau1, 1.0)
    for k in range(499):
        u[k + 1] = a * u[k] + b * i[k]
    assert np.allclose(sim.v1, u, atol=1e-12)
    assert np.allclose(sim.v, ocv(sim.soc) - p.r0 * i - u, atol=1e-12)
    assert p.c1 == pytest.approx(p.tau1 / p.r1)


def test_fit_thevenin_recovers_parameters_noiseless():
    ocv, dt = make_ocv(), 1.0
    truth = TheveninParams(0.045, 0.025, 60.0)
    i = make_current(8000)
    sim = Thevenin1RC(truth, ocv).simulate(i, dt, 0.95, 5.0)
    fit = fit_thevenin_1rc(i, sim.v, sim.soc, ocv, dt)
    assert fit.params.r0 == pytest.approx(truth.r0, rel=1e-3)
    assert fit.params.r1 == pytest.approx(truth.r1, rel=1e-3)
    assert fit.params.tau1 == pytest.approx(truth.tau1, rel=1e-3)
    assert fit.rmse_v < 1e-6
    assert not fit.at_bounds


def test_fit_thevenin_with_noise_and_1rc_beats_rint():
    ocv, dt = make_ocv(), 1.0
    truth = TheveninParams(0.045, 0.025, 60.0)
    i = make_current(8000, seed=3)
    sim = Thevenin1RC(truth, ocv).simulate(i, dt, 0.95, 5.0)
    v = sim.v + np.random.default_rng(2).normal(0, 1e-3, sim.v.size)
    fit = fit_thevenin_1rc(i, v, sim.soc, ocv, dt)
    assert fit.params.tau1 == pytest.approx(60.0, rel=0.15)
    assert fit.params.r0 == pytest.approx(0.045, rel=0.1)
    assert fit.rmse_v == pytest.approx(1e-3, rel=0.15)   # ~ nível de ruído
    rint = fit_rint(i, v, sim.soc, ocv)
    assert fit.rmse_v < rint.rmse_v                      # base do Gate G1


def test_validation_on_different_profile():
    """Identificar num perfil e validar em outro (protocolo de Hu et al.)."""
    ocv, dt = make_ocv(), 1.0
    truth = TheveninParams(0.045, 0.025, 60.0)
    model = Thevenin1RC(truth, ocv)
    i1, i2 = make_current(6000, seed=10), make_current(6000, seed=11)
    s1, s2 = (model.simulate(x, dt, 0.95, 5.0) for x in (i1, i2))
    fit = fit_thevenin_1rc(i1, s1.v, s1.soc, ocv, dt)
    val = evaluate_voltage(Thevenin1RC(fit.params, ocv), i2, s2.v, s2.soc, dt)
    assert val["rmse"] < 1e-5
    assert val["residual"].shape == s2.v.shape


# ------------------------------------------------------------- métricas
def test_metrics_known_values():
    ref = np.array([1.0, 2.0, 3.0, 4.0])
    est = np.array([1.0, 2.0, 3.0, 8.0])
    assert mae(ref, est) == pytest.approx(1.0)
    assert rmse(ref, est) == pytest.approx(2.0)
    assert max_abs_error(ref, est) == pytest.approx(4.0)


def test_inputs_are_validated():
    with pytest.raises(ValueError):
        integrate_soc([1.0, np.nan], 1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        RintModel(RintParams(0.05), make_ocv()).voltage([1, 2], [0.5])


# ----------------------------------------------------------- persistência
def test_params_roundtrip(tmp_path):
    p = TheveninParams(0.04, 0.02, 50.0)
    f = tmp_path / "p.json"
    save_params(str(f), "thevenin_1rc", p, {"perfil": "DST", "T_C": 25})
    q, meta = load_params(str(f))
    assert q == p and meta["perfil"] == "DST"
