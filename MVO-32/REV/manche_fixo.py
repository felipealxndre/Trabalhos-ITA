import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Configuração Estética
sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

# --- 1. DADOS DE ENTRADA (Extraídos das Imagens) ---
# Posições de CG (% CMA)
h_cg1 = 0.2015
h_cg2 = 0.2299


# Dados completos
# # Voo 1 (CG Dianteiro)
# Vi_1 = np.array([180, 170, 160, 150, 140, 140, 150, 160, 170, 178])
# de_1 = np.array([25.0, 25.3, 25.6, 26.1, 26.65, 26.9, 26.4, 25.9, 25.7, 25.4]) * 0.0328084

# # Voo 2 (CG Traseiro)
# Vi_2 = np.array([180, 170, 160, 150, 140, 140, 150, 160, 170, 180])
# de_2 = np.array([25.0, 25.4, 25.7, 26.2, 26.6, 26.5, 26.3, 25.9, 25.5, 25.1]) * 0.0328084


# Dados com média
# Voo 1 (CG Dianteiro)
Vi_1 = np.array([180, 170, 160, 150, 140])
de_1 = np.array([25.00, 25.50, 25.75, 26.25, 26.78]) * 0.0328084

# Voo 2 (CG Traseiro)
Vi_2 = np.array([180, 170, 160, 150, 140])
de_2 = np.array([25.05, 25.45, 25.80, 26.25, 26.55]) * 0.0328084


# Gradientes de estabilidade (d_delta / d_V) [ft/kt]
fit_1 = np.polyfit(Vi_1, de_1, 1)
grad_V1 = fit_1[0] # Slope Voo 1

fit_2 = np.polyfit(Vi_2, de_2, 1)
grad_V2 = fit_2[0] # Slope Voo 2

# Extrapolação para Ponto Neutro
CGs = np.array([h_cg1, h_cg2])
Params = np.array([-grad_V1, -grad_V2])

fit_np = np.polyfit(CGs, Params, 1)
slope_np = fit_np[0]
intercept_np = fit_np[1]

# Ponto Neutro - ponto onde o parâmetro de estabilidade é zero
hn_vel = -intercept_np / slope_np


# Margens Estáticas Fixo
SM_fixo1 = hn_vel - h_cg1
SM_fixo2 = hn_vel - h_cg2

print(f"Ponto Neutro Manche Fixo (hn): {hn_vel*100:.2f} % CMA")
print(f"Margem Estática Fixo Voo 1: {SM_fixo1*100:.2f} %")
print(f"Margem Estática Fixo Voo 2: {SM_fixo2*100:.2f} %")

palette = sns.color_palette("Set2", n_colors=2)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Manche vs Velocidade
ax1.plot(Vi_1, de_1, 'o-', label=f'Voo 1 (CG {h_cg1*100:.2f}%)', color=palette[0])
ax1.plot(Vi_2, de_2, 's-', label=f'Voo 2 (CG {h_cg2*100:.2f}%)', color=palette[1])
v_range = np.linspace(140, 180, 10)
# ax1.plot(v_range, np.polyval(fit_1, v_range), '--', alpha=0.5, color=palette[0])
# ax1.plot(v_range, np.polyval(fit_2, v_range), '--', alpha=0.5, color=palette[1])
ax1.set_xlabel('Velocidade Indicada ($V_{IAS}$) [kt]', fontsize=13)
ax1.set_ylabel('Deflexão do profundor [ft]', fontsize=13)
ax1.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.grid(axis='x', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.legend()
sns.despine(ax=ax1)
for spine in ax1.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

# Plot 2: Extrapolação do Ponto Neutro
ax2.plot(CGs*100, Params, 'ko-', color='r', markersize=8, label='Gradiente Experimental')
limit_viz = 300 if abs(hn_vel*100) > 300 else hn_vel*100 + 10
x_min = min(min(CGs)*100, limit_viz) - 5
x_max = max(max(CGs)*100, limit_viz) + 5

cg_range = np.linspace(15, x_max - 10, 100)
param_range = np.polyval(fit_np, cg_range/100)

ax2.plot(cg_range, param_range, 'r--', label='Extrapolação Linear')
ax2.axhline(0, color='gray', linewidth=0.5)
if abs(hn_vel*100) < 500:
    ax2.plot(hn_vel*100, 0, marker='s', color='g', markersize=8)
ax2.set_xlabel("Posição do CG (% CMA)", fontsize=13)
ax2.set_ylabel(r"Derivada $-d\delta_e / dV$ [ft/kt]", fontsize=13)
ax2.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax2.grid(axis='x', linestyle='--', linewidth=0.8, alpha=0.6)
ax2.legend()
sns.despine(ax=ax2)
for spine in ax2.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

plt.tight_layout()
plt.savefig('analise_manche_fixo.png')
plt.show()