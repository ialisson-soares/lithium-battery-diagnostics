# Proposta de Reestruturação do TCC (v1.2)

---

## 1. Papel de PF1

PF1 deve responder essencialmente à seguinte pergunta:

> **É tecnicamente viável modelar a bateria e estimar seus estados, em especial o SoC, utilizando uma abordagem suficientemente simples para posterior implementação em uma PoC embarcada?**

PF1 passa a ser **software + dados + especificação da PoC**. Não deve haver compromisso com a construção de hardware neste momento.

### Cadeia Mínima do Escopo
$$\text{Fundamentação} \longrightarrow \text{Dataset} \longrightarrow \text{Modelo} \longrightarrow \text{Estimador} \longrightarrow \text{Avaliação} \longrightarrow \text{Especificação da PoC}$$

---

### Papel de PF2
PF2 responde a outra pergunta fundamental:

> **O método selecionado em PF1 pode ser implementado e avaliado em uma plataforma física de baixo custo?**

---

## 2. O que deve permanecer obrigatoriamente em PF1

O escopo de PF1 fica restrito a apenas **seis blocos essenciais**:

### A. Problema da bateria e fundamentação
Revisão dirigida, sem tentar produzir um estado da arte extenso:
- Comportamento da bateria Li-ion;
- Tensão ($V$), corrente ($I$), temperatura ($T$);
- Estado de Carga ($SoC$) e Estado de Saúde ($SoH$);
- Tensão em Circuito Aberto vs. SoC ($OCV–SoC$);
- Capacidade e resistência interna;
- Processos de degradação;
- Função de um BMS e a diferença entre monitoramento passivo e BMS completo.

*O resultado deve ser suficiente para justificar as decisões posteriores, não uma revisão enciclopédica.*

### B. Um modelo de bateria
Adoção direta do modelo:
$$\boxed{\text{ECM 1RC}}$$

* $R_\text{int}$ servirá apenas como explicação introdutória/baseline.
* Não haverá comparação sistemática temporal entre $R_\text{int}$, $1RC$, $2RC$ e modelos eletroquímicos. O $1RC$ é suficientemente simples para implementação posterior e representativo para o objetivo acadêmico do trabalho.

### C. Estado de Carga (SoC)
Manutenção de dois métodos para comparação técnica significativa:
1. $\boxed{\text{Coulomb Counting}}$ *(baseline)*
2. $\boxed{\text{ECM 1RC + EKF}}$ *(método principal)*

### D. Estado de Saúde (SoH) — Escopo Reduzido
Sem estimadores sofisticados de $SoH$ em PF1. Utilização da relação direta de capacidade:

$$SoH_Q = \frac{Q_\text{atual}}{Q_\text{referência}}$$

utilizando dados públicos de envelhecimento para demonstrar a compreensão de perda de capacidade, relação com envelhecimento e extração do indicador a partir de um dataset. A estimação online de $SoH$ fica para PF2 ou extensões.

### E. Dados sintéticos + Open Dataset
Papéis distintos para cada fonte:
* **Dados Sintéticos:** Não serão usados para validar o modelo que os gerou. Servirão para verificar comportamentos controlados (erro inicial de $SoC$, bias de corrente, ruído de tensão, mudança de capacidade).
* **Open Datasets (Evidência Externa):**
  * **CALCE:** Dataset principal para a validação de $SoC$.
  * **NASA Aging Dataset:** Subconjunto utilizado apenas para demonstração preliminar do $SoH_Q$.

### F. Especificação da PoC
PF1 é concluído com um documento detalhado para que PF2 inicie imediatamente:
- Escolha da bateria/célula e microcontrolador;
- Sensores e circuito de aquisição de $V, I, T$;
- Frequência de amostragem e taxas de interface;
- Estratégias de armazenamento, telemetria e alimentação;
- Proteção da bancada experimental;
- Algoritmo que será portado e protocolo experimental.

---

## 3. Cronograma Realista para PF1

| Período | Atividade | Produto |
| :--- | :--- | :--- |
| **22–27/out** | Revisão dirigida + delimitação | Definição do problema, requisitos e referências |
| **28/out–03/nov** | Seleção e preparação do CALCE | Notebook de importação e exploração |
| **04–10/nov** | ECM 1RC + $OCV–SoC$ | Modelo executável |
| **11–17/nov** | Coulomb Counting + EKF | Comparação de $SoC$ |
| **18–23/nov** | Perturbações sintéticas + $SoH$ preliminar NASA | Resultados experimentais |
| **24–30/nov** | Consolidação + especificação da PoC | Relatório PF1 + plano para PF2 |

---

## 4. Entregáveis Mínimos de PF1

1. **PF1-D1 — Fundamentação e Trabalhos Relacionados:** Documento sintético cobrindo bateria Li-ion, $SoC$, $SoH$, ECM, Coulomb Counting, EKF e trabalhos correlatos essenciais.
2. **PF1-D2 — Dataset e Pipeline:** Notebooks/scripts para carregamento do dataset CALCE, identificação de $V, I, T, t$, pré-processamento, geração de gráficos básicos e rastreabilidade dos arquivos.
3. **PF1-D3 — Modelo:** Implementação do $ECM\ 1RC$ acompanhada da documentação e parametrização.
4. **PF1-D4 — Estimação de SoC:** Implementação e comparação entre Coulomb Counting e EKF, utilizando as métricas mínimas $MAE_{SoC}$ e $RMSE_{SoC}$.
5. **PF1-D5 — SoH e Análise de Sensibilidade:** Demonstração do $SoH$ baseada em capacidade (subconjunto NASA) e análise de sensibilidade a erros iniciais de $SoC$, bias de corrente e ruído de medição.
6. **PF1-D6 — Especificação de PF2:** Documentação da arquitetura da PoC e plano experimental detalhado.

---

## 5. Resultado Esperado de PF1

Ao final do período, o trabalho entregará a seguinte consolidação:

> *"Foi implementado um modelo equivalente 1RC de bateria Li-ion. A estimação de SoC por contagem de Coulomb foi utilizada como baseline e comparada a um estimador EKF baseado no modelo. Os métodos foram inicialmente submetidos a perturbações controladas e posteriormente avaliados utilizando dados experimentais públicos. Foi realizada ainda uma análise preliminar de SoH baseada em capacidade utilizando dados públicos de envelhecimento. Com base nesses resultados, foi especificada a plataforma experimental a ser implementada em PF2."*

---

## 6. Itens Retirados Explicitamente de PF1

| Item | Decisão |
| :--- | :--- |
| Digital Twin | Retirar |
| RUL (Remaining Useful Life) | Retirar |
| EIS (Espectroscopia de Impedância) | Retirar |
| Machine Learning | Retirar |
| Modelo 2RC | Não implementar *(salvo necessidade evidente)* |
| Comparação extensa de modelos | Retirar |
| Múltiplos open datasets | Retirar |
| Hardware físico / PCB | Escopo de PF2 |
| Firmware embarcado | Escopo de PF2 |
| Telemetria | Escopo de PF2 |
| Caracterização de sensores | Escopo de PF2 |
| $SoH$ online | Escopo de PF2 / Extensão |
| Diagnóstico avançado | Escopo de PF2 / Extensão |

---

## 7. Estrutura de PF2 (Foco na PoC)

PF2 fica organizado em cinco blocos principais:

* **PF2-A — Projeto da PoC:** Implementação física da arquitetura especificada:
  $$\text{Bateria} \longrightarrow \text{Sensores } (V/I/T) \longrightarrow \text{MCU} \longrightarrow \text{Estimação} \longrightarrow \text{Telemetria}$$
* **PF2-B — Caracterização da Instrumentação:** Levantamento prévio obrigatório de erro de tensão, erro de corrente, erro de temperatura, ruído, offset e frequência efetiva de amostragem.
* **PF2-C — Port do Algoritmo:** Implementação no MCU na seguinte ordem sequencial:
  1. Aquisição de dados;
  2. Coulomb Counting;
  3. Modelo ECM;
  4. EKF *(se os resultados de PF1 justificarem)*.
* **PF2-D — Experimentação Física:** Produção de dataset próprio contendo $t, V, I, T, \widehat{SoC}$ através de ciclos seguros de carga e descarga.
* **PF2-E — Comparação Final:** Estrutura completa de validação:
  $$\boxed{\text{Sintético}} \longrightarrow \boxed{\text{Open Dataset}} \longrightarrow \boxed{\text{PoC Física}}$$

---

## 8. Cronograma Conceitual de PF2

1. Aquisição dos componentes e montagem
2. Firmware de aquisição
3. Calibração e caracterização
4. Coulomb Counting embarcado
5. Modelo / EKF embarcado
6. Telemetria
7. Experimentos de bancada
8. Comparação com os resultados de PF1
9. Análise de limitações
10. Redação final do TCC

---

## 9. Expectativa para a Estimação de SoH

A validação de $SoC$ e $SoH$ assume abordagens distintas de acordo com a viabilidade experimental de tempo:

$$\text{SoC:} \quad \text{Sintético} + \text{Open Dataset} + \text{PoC Física}$$

$$\text{SoH:} \quad \text{Open Dataset} + \text{Análise Offline } (+ \text{Indicador preliminar na PoC, se viável})$$

---

## Revisão da Estratégia Geral

* **PF1:** *Modelagem e avaliação de métodos de estimação de estado de bateria.*
* **PF2:** *Implementação e validação experimental em uma plataforma embarcada conectada.*