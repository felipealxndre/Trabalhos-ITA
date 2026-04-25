# AH-1S Mission Performance — PRJ-91 / Laboratório nº 3

Simulação computacional de desempenho e consumo de combustível do helicóptero **Bell AH-1S Cobra** em uma missão de transporte. O código implementa o modelo aerodinâmico do voo pairado e do voo à frente baseado na **Teoria da Quantidade de Movimento (QDM/Froude)** combinada com a **Teoria do Elemento de Pá (TEP)**, conforme apresentado nas referências (Leishman, Prouty, Johnson).

---

## Sumário

1. [Estrutura do projeto](#1-estrutura-do-projeto)
2. [Como rodar](#2-como-rodar)
3. [Cenários simulados](#3-cenários-simulados)
4. [Modelo aerodinâmico](#4-modelo-aerodinâmico)
5. [Metodologia da simulação](#5-metodologia-da-simulação)
6. [Hipóteses adotadas](#6-hipóteses-adotadas)
7. [Justificativa das decisões de implementação](#7-justificativa-das-decisões-de-implementação)
8. [Saídas](#8-saídas)
9. [Limitações conhecidas](#9-limitações-conhecidas)
10. [Referências](#10-referências)

---

## 1. Estrutura do projeto

```
LAB-02/
├── main.py            # ponto de entrada — roda os 4 cenários
├── config.py          # dados do AH-1S, parâmetros de missão, ISA, conversões
├── aero.py            # inflow de Glauert, potência, efeito solo, curva PN(V0), velocidades
├── simulation.py      # fases da missão e simulação completa
├── report.py          # tabelas (console + CSV), balanço de potência, polar de velocidades
├── README.md          # este documento
├── requirements.txt   # dependências Python
└── output/
    ├── results/       # CSVs gerados (criado no runtime)
    └── plots/         # PNGs gerados (criado no runtime)
```

### Princípios de organização

- **5 módulos coesos**: o código foi consolidado de 16 arquivos dispersos para 5 módulos com responsabilidades claras — configuração, aerodinâmica, simulação, relatório e ponto de entrada.
- **Imutabilidade**: os dados de entrada (`Aircraft`, `Rotor`, `Mission`) são `frozen=True` dataclasses para evitar mutações acidentais.
- **Funções puras**: as rotinas de cálculo aerodinâmico são puras (mesmo input → mesmo output, sem efeitos colaterais), facilitando testes e validação.
- **Tudo em SI internamente**: a entrada/saída lida com unidades inglesas (lb, ft, NM, kt) por imposição do enunciado, mas todo o processamento ocorre em SI. As conversões ficam confinadas em `config.py`.

---

## 2. Como rodar

### Dependências

```bash
pip install -r requirements.txt
```

Bibliotecas usadas: `numpy`, `scipy`, `matplotlib`, `pandas`.

### Execução

```bash
python main.py
```

O script roda os 4 cenários sequencialmente. Para cada cenário:

1. Imprime no console a tabela de consumo por fase (kW por parcela, lb por fase).
2. Imprime a tabela de velocidades características (Vy, VAM, VDM, VMR, VH).
3. Imprime verificação de potência (PN vs PU por fase).
4. Salva CSVs em `output/results/`.
5. Gera plots de balanço de potência e polar de velocidades para fases 2–5 em `output/plots/`.

---

## 3. Cenários simulados

Todos os cenários partem da **massa máxima de decolagem (10.000 lb)**, atmosfera **ISA+20**, e seguem a sequência fixa de 6 fases:

| # | Fase | Duração / Distância | Característica chave |
|---|---|---|---|
| 1 | Pairado IGE inicial | 5 min, h = 6 ft, Zp = 0 | Efeito solo reduz Vi |
| 2 | Subida em Vy | 0 ft → 5.000 ft | Vz > 0, Cp inclui λc·CT, correção de λc no inflow |
| 3 | Cruzeiro nivelado | em Zp = 5.000 ft | Voa em **VDM** |
| 4 | Reserva | 30 min em Zp = 5.000 ft | Voa em **VAM** |
| 5 | Descida em Vy | 5.000 ft → 0 ft | Vz < 0, Pi calculada como nivelado (λc=0 no Glauert) |
| 6 | Pairado IGE final | 5 min, h = 6 ft, Zp = 0 | Efeito solo |

### Os quatro cenários

| # | Cenário | Distância | Vz subida | Vento de proa |
|---|---|---|---|---|
| 1 | Baseline | 400 NM | 1.000 ft/min | 0 kt |
| 2 | Vento de proa | 400 NM | 1.000 ft/min | **15 kt** |
| 3 | Subida rápida | 400 NM | **2.000 ft/min** | 0 kt |
| 4 | Distância maior | **440 NM** | 1.000 ft/min | 0 kt |

---

## 4. Modelo aerodinâmico

### 4.1 Coeficiente de potência total

A equação central, conforme o enunciado do laboratório, é:

```
Cp = ki · CT² / [2·√(μ² + (λc + λi)²)]    (induzida)
   + (σ·Cd0/8) · (1 + 4,65·μ²)             (perfil)
   + (1/2) · (f/A) · μ³                    (parasita)
   + Cp_misc                               (miscelânea)
   + λc · CT                               (subida/descida)
```

Onde:

| Símbolo | Significado | Origem |
|---|---|---|
| `CT` | Coeficiente de tração | `T / [ρ·A·(ΩR)²]` |
| `μ` | Razão de avanço | `V₀ / (ΩR)` |
| `λi` | Inflow induzido | Glauert (implícita) |
| `λc` | Inflow de subida | `Vz / (ΩR)` |
| `σ` | Solidez do rotor | `b·c / (π·R)` |
| `Cd0` | Arrasto de perfil | Constante |
| `f` | Placa plana equivalente | Catálogo Prouty |
| `ki` | Correção de perdas QDM | 1,15 |
| `ηm` | Eficiência mecânica | 0,85 |

### 4.2 Inflow induzido (Glauert)

A razão $\lambda_i$ é raiz de uma equação implícita:

$$\lambda_i = \frac{C_T}{2\sqrt{\mu^2 + (\lambda_c + \lambda_i)^2}}$$

Resolvida em `aero.py` por **iteração de ponto fixo com relaxação**:

$$\lambda_i^{(n+1)} = \lambda_i^{(n)} + \omega \cdot [f(\lambda_i^{(n)}) - \lambda_i^{(n)}]$$

com:
- Chute inicial: $\lambda_i^{(0)} = \sqrt{C_T/2}$ (valor exato do pairado).
- Fator de relaxação: $\omega = 0{,}5$ (estabiliza a convergência em $\mu$ pequeno).
- Tolerância: $10^{-8}$.
- Convergência típica: 10-30 iterações.

A relaxação é necessária porque, próximo ao pairado ($\mu \to 0$), a iteração simples pode oscilar em torno da solução. Com $\omega = 0{,}5$, garante-se convergência monotônica em todo o envelope.

### 4.3 Correção de λc na subida

Na fase de subida, a curva de potência PN(V₀) é calculada **incluindo o efeito de λc no cálculo iterativo do inflow** (Eq. de Glauert completa). Isso reduz a potência induzida, especialmente em baixas velocidades, deslocando o ponto de potência mínima (VAM) para a esquerda e alterando ligeiramente Vy. O efeito é mais pronunciado com Vz maiores — no Caso 3 (2.000 ft/min), Vy cai de ~72,7 para ~70,1 kt.

Nos gráficos de balanço de potência da subida, a curva PN de voo nivelado (sem correção λc) é sobreposta para visualizar o efeito da correção.

### 4.4 Tratamento da descida

Na descida, conforme instrução do laboratório:
- A **potência induzida é calculada como em voo nivelado** (λc=0 no cálculo do inflow de Glauert), desprezando a influência da velocidade vertical sobre Vi.
- A **parcela λc·CT é mantida** com o Vz real (negativo), contabilizando o ganho de energia potencial que reduz a potência total requerida.

### 4.5 Potência de miscelânea

A miscelânea engloba o rotor de cauda, acessórios e perdas de transmissão. Calculada a partir da eficiência mecânica:

$$P_{misc} = \left(\frac{1}{\eta_m} - 1\right) \cdot P_{eR}$$

onde $P_{eR}$ é a potência no eixo do rotor principal (induzida + perfil + parasita). Com $\eta_m = 0{,}85$, isso resulta em $C_{P,misc} \approx 0{,}176 \cdot C_{P,eR}$.

### 4.6 Efeito solo

A redução de potência induzida no pairado IGE é calculada via **ajuste empírico** aos pontos da curva tracejada da Figura 3 do Prouty (full-scale flight test results, p. 66). O ajuste usa a função:

$$\frac{V_{1,IGE}}{V_{1,OGE}} = 1 - a \cdot e^{-b \cdot (z/D)}$$

Os parâmetros $a$ e $b$ são determinados por mínimos quadrados não-lineares (`scipy.optimize.curve_fit`) a partir de 11 pontos lidos do gráfico, cobrindo $z/D$ de 0,10 a 1,20.

No código, implementado como um **fator multiplicativo** sobre a parcela induzida do $C_p$.

### 4.7 Atmosfera ISA + ΔT

Implementação em `config.py`:

- **Temperatura ISA**: $T_{ISA}(Z_p) = 288{,}15 - 0{,}0065 \cdot Z_p[m]$ K
- **Temperatura efetiva**: $T = T_{ISA} + \Delta T$
- **Pressão**: $p = p_0 \cdot (T_{ISA}/T_0)^{g/(R \cdot L)}$ (gradiente padrão)
- **Densidade**: $\rho = p / (R \cdot T)$ (gás ideal com T efetiva)

### 4.8 Potência disponível com efeito Ram

$$P_U(V_0) = P_{U,0} \cdot \left(1 + k_{ram} \cdot \frac{\rho V_0^2}{2 P_{atm}}\right)$$

Adotou-se $k_{ram} = 0{,}7$, valor típico para admissões de helicópteros sem otimização agressiva. Sem efeito Ram, $V_y$ coincidiria exatamente com $V_{AM}$, o que seria uma simplificação imprecisa.

---

## 5. Metodologia da simulação

### 5.1 Definição das velocidades características

Implementadas em `aero.py`:

| Velocidade | Definição matemática | Significado físico |
|---|---|---|
| $V_{AM}$ | $\arg\min_{V_0} P_N(V_0)$ | Mínimo da curva PN — máxima autonomia |
| $V_{DM}$ | $\arg\min_{V_0} [P_N(V_0) / (V_0 - V_w)]$ | Tangente a partir de Vw — máximo alcance |
| $V_y$ | $\arg\max_{V_0} [P_U(V_0) - P_{N,aero}(V_0)]$ | Máximo excesso de potência — máxima razão de subida |
| $V_{MR}$ | $= V_{AM}$ | Menor razão de descida em autorrotação |
| $V_H$ | $\max\{V_0 : P_N(V_0) \leq P_U(V_0)\}$ | Velocidade máxima nivelada |
| $V_{z,max}$ | $\frac{(P_U - P_{N,aero}) \cdot \eta_m}{W}$ em $V_y$ | Razão de subida atingível |

**Observações importantes:**

- **VDM com vento**: a velocidade de máximo alcance é calculada minimizando $P_N / V_{ground}$ onde $V_{ground} = V_{TAS} - V_w$. Graficamente, corresponde à reta tangente que parte de $V_w$ no eixo das abscissas (não da origem).

- **Vy na subida**: quando a curva é computada com a correção de λc (Vz no inflow), a potência induzida diminui em baixas velocidades, e $P_{N,aero}$ (sem a parcela de subida constante) é usada para extrair Vy. Como a parcela de subida $\lambda_c \cdot C_T$ é constante com V₀, ela não afeta o argmax.

- A discretização de V₀ é feita em **2.000 pontos** entre 0,5 m/s e 100 m/s (~194 kt), com precisão ~0,1 kt no resultado final.

### 5.2 Cálculo de tempo e distância por fase

| Fase | Tempo | Distância |
|---|---|---|
| Pairado | Fixo (5 min) | 0 |
| Subida | $\Delta h / V_z$ | $V_y \cdot t$ |
| Cruzeiro | $d / V_{ground}$ | Fixa (400 ou 440 NM) |
| Reserva | Fixo (30 min) | $V_{AM} \cdot t$ |
| Descida | $\Delta h / |V_z|$ | $V_y \cdot t$ |

No cruzeiro com vento de proa:

$$V_{ground} = V_{TAS} - V_{vento} \qquad t = d_{NM} / V_{ground}$$

### 5.3 Cálculo de combustível por fase

Consumo específico convertido para SI:

$$\text{SFC}_{SI} = 0{,}458 \, \frac{\text{lb}}{\text{hp·hr}} \cdot \frac{0{,}45359 \, \text{kg/lb}}{745{,}7 \, \text{W/hp} \cdot 3600 \, \text{s/hr}} = 7{,}735 \times 10^{-8} \, \frac{\text{kg}}{\text{W·s}}$$

Para cada fase:

$$m_{comb} = \text{SFC}_{SI} \cdot P_{N,total} \cdot t$$

### 5.4 Iteração de peso médio

Em fases longas (cruzeiro, reserva), a massa varia significativamente. Tratamento:

1. Estimativa inicial: $P_{N}$ calculada com massa inicial → consumo estimado.
2. Massa média: $\bar{m} = m_{inicial} - 0{,}5 \cdot m_{comb,estimado}$.
3. Recalcula $P_N$ com $\bar{m}$ → novo consumo.
4. Itera 3 vezes (suficiente para convergência < 0,1%).

Implementação no helper `_avg_weight_iter()` em `simulation.py`.

### 5.5 Sequência da simulação

```
mass = max_takeoff
para cada fase em [pairado_inicial, subida, cruzeiro, reserva, descida, pairado_final]:
    1. Calcula condições atmosféricas (Zp, ISA+ΔT) → ρ
    2. Gera curva PN(V0) com correção apropriada:
       - Subida: inclui λc no inflow (curva corrigida)
       - Descida: λc=0 no inflow (como nivelado)
       - Demais: Vz=0
    3. Extrai velocidade característica relevante (Vy, VDM ou VAM)
    4. Aplica iteração de peso médio (se relevante)
    5. Calcula PN total e suas parcelas
    6. Calcula consumo de combustível
    7. Atualiza massa: mass -= fuel_burned
    8. Registra resultados
```

---

## 6. Hipóteses adotadas

### 6.1 Hipóteses do modelo aerodinâmico

**Da QDM (Teoria da Quantidade de Movimento):**
1. Fluido incompressível (válido para $M < 0{,}3$).
2. Escoamento estacionário e sem viscosidade (perfil sem arrasto — corrigido pela TEP).
3. Velocidade induzida uniforme ao longo do disco.
4. Hipótese de Glauert para o voo à frente.
5. Ângulo de ataque do disco $\alpha_R \approx 0$.

**Da TEP:**
6. $C_{d_0}$ constante ao longo da envergadura e do azimute.
7. Pequenos ângulos: $\sin x \approx x$, $\cos x \approx 1$.
8. Pás retangulares com corda constante.
9. Sem estol da pá que recua nem efeitos de compressibilidade (válido até $\mu \approx 0{,}3$).

**Da fórmula CPp = (σCd₀/8)(1 + 4,65μ²):**
10. Inclui contribuições do arrasto de rotação, translação azimutal e radial.
11. Coeficiente 4,65 é uma aproximação consagrada.

### 6.2 Hipóteses do modelo de motor

12. **Potência disponível** segue o modelo Ram: $P_U(V_0) = P_{U,0} \cdot (1 + 0{,}7 \cdot \rho V_0^2/(2P_{atm}))$.
13. **Consumo específico SFC = 0,458 lb/(hp·hr)** constante.
14. **Eficiência mecânica $\eta_m = 0{,}85$** constante em todas as condições.

### 6.3 Hipóteses operacionais

15. **Vento constante** em rota (não variando com altitude nem ao longo da trajetória).
16. **Pairado IGE** com $h = 6$ ft fixo; **z/D = 6/44 = 0,136**.
17. **Subida e descida em Vy** com razão constante (1.000 ou 2.000 ft/min).
18. **Cruzeiro em VDM**, reserva em VAM.
19. **Densidade média da Zp média** para subida e descida (Zp_média = 2.500 ft).
20. **Variação de Vi com Vz desprezada na descida** (lab pede explicitamente). Mantida na subida, com efeito de λc no cálculo iterativo do inflow.
21. **VMR = VAM** no modelo simplificado (sem modelagem de autorrotação completa).

---

## 7. Justificativa das decisões de implementação

### 7.1 Por que consolidar em 5 módulos?

O projeto original continha 16 arquivos Python distribuídos em 6 diretórios. A dispersão excessiva dificultava a leitura e o acompanhamento do fluxo de dados. A consolidação em 5 módulos (`config.py`, `aero.py`, `simulation.py`, `report.py`, `main.py`) mantém a separação de responsabilidades mas reduz drasticamente a fragmentação, facilitando revisão e manutenção.

### 7.2 Por que iteração de ponto fixo com relaxação para λi?

A equação de Glauert é implícita em $\lambda_i$. O método de ponto fixo com relaxação ($\omega = 0{,}5$) é robusto em todo o envelope, simples e sem dependência de derivadas. Convergência típica em 10-30 iterações.

### 7.3 Por que ajuste empírico ao Prouty em vez de Cheeseman-Bennett?

A fórmula clássica de Cheeseman-Bennett só vale para $z/R > 0{,}5$ ($z/D > 0{,}25$). Para $h = 6$ ft / $R = 22$ ft, temos $z/R = 0{,}273$, **fora da faixa de validade**. O ajuste à curva da Fig. 3 do Prouty, baseada em ensaios em voo full-scale, é mais fiel.

### 7.4 Por que correção de λc na curva de subida?

Conforme o relatório de referência do laboratório (Lab03_PRJ91, 2025, Seção 2.4), na subida o aumento do fluxo mássico através do rotor reduz a potência induzida, especialmente em baixas velocidades. Incluir λc no cálculo iterativo do inflow para a curva PN(V₀) da fase de subida reflete este efeito e desloca a VAM para a esquerda. Nos gráficos, a curva nivelada (sem correção) é sobreposta para comparação.

### 7.5 Por que VDM com tangente a partir de Vw?

O alcance máximo depende da velocidade em relação ao **solo**, não ao ar. Com vento de proa de magnitude $V_w$, a tangente à curva PN parte do ponto $(V_w, 0)$ no eixo das abscissas, resultando em VDM maior que no caso sem vento. Isso é implementado minimizando $P_N / (V_{TAS} - V_w)$.

### 7.6 Por que polar de velocidades?

A polar de velocidades Vz(V₀) mostra graficamente a capacidade de subida/descida da aeronave em cada velocidade. A partir dela identificam-se VvM (= Vy, razão de subida máxima), Vvm (menor razão de descida), e VH (onde Vz cruza zero). Este gráfico é gerado para cada fase (2–5) e complementa o balanço de potência na documentação dos resultados.

### 7.7 Por que peso médio em vez de integração numérica?

A massa varia ~15% no cruzeiro de 400 NM. Peso médio com 3 iterações tem erro < 0,1% comparado à integração numérica, com ~10× menos esforço computacional.

### 7.8 Por que SI internamente?

O enunciado mistura sistemas de unidades. Tratar tudo em SI elimina ambiguidades, centraliza conversões em `config.py`, e permite validação por análise dimensional direta.

---

## 8. Saídas

### 8.1 Console

Para cada cenário, são impressos:

- **Tabela de consumo por fase**: duração, distância, massa inicial e final, combustível gasto, velocidade utilizada, e parcelas de potência (induzida, perfil, parasita, miscelânea, subida/descida, total) em kW.
- **Tabela de velocidades características** (fases 2 a 5): Vy, VAM, VDM, VMR, VH em kt e Rz_max em ft/min.
- **Resumo**: massa inicial, massa final, combustível total, capacidade do tanque, margem.
- **Verificação de potência**: compara PN vs PU para cada fase, indicando [OK] ou [EXCEDE PU].

### 8.2 Arquivos CSV

Em `output/results/`:

- `Caso_X_consumo.csv`: tabela completa de consumo por fase.
- `Caso_X_velocidades.csv`: velocidades características.

### 8.3 Plots

Em `output/plots/`, para cada cenário são gerados **8 gráficos** (2 por fase, fases 2–5):

**Balanço de potência** (`*_balanco.png`):
- Curva PN total (preta sólida) com parcelas individuais (tracejadas).
- Linha de PU com efeito Ram (vermelha pontilhada).
- Marcadores em Vy (laranja), VAM (verde), VDM (azul), VH (roxo).
- Reta tangente origem-VDM (azul pontilhada; parte de Vw se houver vento).
- Na subida: componente de subida (marrom) e curva nivelada (cinza) para comparação.
- Na descida: componente de descida negativa (marrom).

**Polar de velocidades** (`*_polar.png`):
- Curva Vz(V₀) = (PU - PN_aero) · ηm / W [ft/min].
- Marcadores em VvM (= Vy), Vvm (mínima razão de descida), VH.
- Na subida: curva com correção de λc no inflow.

Total: **32 gráficos** (8 por cenário × 4 cenários).

---

## 9. Limitações conhecidas

### 9.1 Limites de validade do modelo aerodinâmico

- **Velocidade máxima**: o termo $(1 + 4{,}65\mu^2)$ para potência de perfil tem precisão de ~5% até $\mu = 0{,}5$. Acima disso, estol da pá que recua e compressibilidade dominam.
- **Pairado e baixa velocidade**: a aproximação de Glauert ainda vale, mas o inflow real é mais não-uniforme. O fator $k_i = 1{,}15$ tenta capturar essas perdas globalmente.
- **Anéis de vórtice**: a QDM não vale entre $V_z/V_{i_0} \in (-2, -0{,}25)$. O código não detecta essa região.

### 9.2 Limitações do modelo de motor

- SFC constante: motores reais têm SFC piorando em regimes parciais e em altitude.
- Modelo Ram simplificado: aproximação para baixos Mach.
- Sem perdas com altitude: a queda de PU é capturada apenas via densidade no modelo Ram.

### 9.3 Precisão numérica

- A varredura PN(V₀) usa **2.000 pontos**, dando precisão de ~0,1 kt para velocidades características.
- A iteração de inflow tem tolerância $10^{-8}$.
- O ajuste do efeito solo tem R² > 0,99.

---

## 10. Referências

1. Leishman, J. Gordon. *Principles of Helicopter Aerodynamics*, 2nd ed. Cambridge University Press, 2006.
2. Prouty, Raymond W. *Helicopter Performance, Stability and Control*. Krieger, 2003.
3. Johnson, W. *Helicopter Theory*. Dover Publications, 1994.
4. Notas de aula PRJ-91 — Fundamentos de Engenharia de Helicópteros, ITA.
5. Lab03_PRJ91 — Relatório do laboratório do ano anterior (2025), referência para metodologia.

---

## Apêndice A — Equações resumo

**Adimensionalização:**
$$C_T = \frac{T}{\rho A (\Omega R)^2} \qquad C_P = \frac{P}{\rho A (\Omega R)^3}$$
$$\mu = \frac{V_0}{\Omega R} \qquad \lambda_i = \frac{V_i}{\Omega R} \qquad \lambda_c = \frac{V_z}{\Omega R}$$

**Pairado (QDM):**
$$V_{i_0} = \sqrt{\frac{mg}{2\rho A}} \qquad P_{i_0} = T \cdot V_{i_0} \qquad \lambda_{i_0} = \sqrt{\frac{C_T}{2}}$$

**Glauert (voo à frente):**
$$\lambda_i = \frac{C_T}{2\sqrt{\mu^2 + (\lambda_c + \lambda_i)^2}}$$

**Coeficiente de potência completo (lab):**
$$C_P = k_i \frac{C_T^2}{2\sqrt{\mu^2 + (\lambda_c + \lambda_i)^2}} + \frac{\sigma C_{d_0}}{8}(1 + 4{,}65\mu^2) + \frac{1}{2}\frac{f}{A}\mu^3 + C_{P,misc} + \lambda_c C_T$$

**Combustível:**
$$m_{comb} = \text{SFC} \cdot P_{N,total} \cdot t$$

**Razão de subida atingível:**
$$V_{z,max} = \frac{(P_U - P_{N,aero}) \cdot \eta_m}{W} \quad \text{em } V_0 = V_y$$
