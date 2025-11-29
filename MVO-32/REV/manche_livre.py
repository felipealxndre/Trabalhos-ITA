import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

# Posições de CG (% CMA)
h_cg1 = 0.2015
h_cg2 = 0.2299


# # Voo 1 (CG Dianteiro)
# Vi_1 = np.array([180, 190, 200, 210, 220, 210, 200, 190, 180])
# Fs_acc1 = np.array([0, -4, -9, -14, -23, -17, -14, -7, 0])

# # Voo 2 (CG Traseiro)
# Vi_2 = np.array([180, 190, 200, 210, 220, 230, 220, 210, 200])
# Fs_acc2 = np.array([0, -5, -9, -15, -24, -31, -19, -14, -10])

# Dados em formato tabular para referência
Vi = np.array([180, 190, 200, 210, 220])
Vi_1 = Vi
Vi_2 = Vi
Fs_acc1 = np.array([0.0, -5.5, -11.5, -17, -23.0])  # np.nan para valores ausentes
Fs_acc2 = np.array([0.0, -5.0, -9.5, -14.5, -19])


# Ajuste linear
fit_livre1 = np.polyfit(Vi_1, Fs_acc1, 1)
grad_P1 = fit_livre1[0] # kgf/kt

fit_livre2 = np.polyfit(Vi_2, Fs_acc2, 1)
grad_P2 = fit_livre2[0] # kgf/kt


CGs = np.array([h_cg1, h_cg2])
Params = np.array([-grad_P1, -grad_P2])
fit_np = np.polyfit(CGs, Params, 1)
slope_np = fit_np[0]
intercept_np = fit_np[1]

hn_livre = -intercept_np / slope_np

# Margens Estáticas Livres
SM_livre1 = hn_livre - h_cg1
SM_livre2 = hn_livre - h_cg2

print(f"Ponto Neutro Manche Livre (hn'): {hn_livre*100:.2f} % CMA")
print(f"Margem Estática Livre Voo 1: {SM_livre1*100:.2f} %")
print(f"Margem Estática Livre Voo 2: {SM_livre2*100:.2f} %")

palette = sns.color_palette("Set2", n_colors=2)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Força no Manche vs Velocidade (invertido)
ax1.plot(Vi_1, -Fs_acc1, 'o-', label=f'Voo 1 (CG {h_cg1*100:.2f}%)', color=palette[0])
ax1.plot(Vi_2, -Fs_acc2, 's-', label=f'Voo 2 (CG {h_cg2*100:.2f}%)', color=palette[1])
v_range = np.linspace(180, 220, 10)
# ax1.plot(v_range, -np.polyval(fit_livre1, v_range), '--', alpha=0.5, color=palette[0])
# ax1.plot(v_range, -np.polyval(fit_livre2, v_range), '--', alpha=0.5, color=palette[1])
ax1.set_xlabel('Velocidade Indicada ($V_{IAS}$) [kt]', fontsize=13)
ax1.set_ylabel('Força no Manche $-F_s$ [kgf]', fontsize=13)
ax1.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.grid(axis='x', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.legend()
sns.despine(ax=ax1)
for spine in ax1.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

# Plot 2: Extrapolação do Ponto Neutro
ax2.plot(CGs*100, Params, 'ko-', color='r', markersize=8, label='Gradiente Experimental')
limit_viz = 300 if abs(hn_livre*100) > 300 else hn_livre*100 + 10
x_min = min(min(CGs)*100, limit_viz) - 5
x_max = max(max(CGs)*100, limit_viz) + 5

cg_range = np.linspace(15, x_max - 10, 100)
param_range = np.polyval(fit_np, cg_range/100)

ax2.plot(cg_range, param_range, 'r--', label='Extrapolação Linear')
ax2.axhline(0, color='gray', linewidth=0.5)
if abs(hn_livre*100) < 500:
    ax2.plot(hn_livre*100, 0, marker='s', color='g', markersize=8)
ax2.set_xlabel("Posição do CG (% CMA)", fontsize=13)
ax2.set_ylabel(r"Derivada $-dF_s / dV$ [kgf/kt]", fontsize=13)
ax2.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax2.grid(axis='x', linestyle='--', linewidth=0.8, alpha=0.6)
ax2.legend()
sns.despine(ax=ax2)
for spine in ax2.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

plt.tight_layout()
plt.savefig('analise_manche_livre.png')
plt.show()