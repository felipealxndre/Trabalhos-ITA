import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Posições de CG (% CMA)
h_cg1 = 0.2015
h_cg2 = 0.2299

# DADOS MANCHE FIXO (Pontos Estabilizados)
Vi_stab1 = np.array([130, 140, 150, 160, 170, 180])
de_stab1 = np.array([27.05, 26.65, 26.10, 25.60, 25.30, 25.00])

Vi_stab2 = np.array([130, 140, 150, 160, 170, 180])
de_stab2 = np.array([27.10, 26.60, 26.20, 25.70, 25.40, 25.00])

# DADOS MANCHE LIVRE (Aceleração - Método Extra 1)
Vi_acc1 = np.array([180, 190, 200, 210, 220])
Fs_acc1 = np.array([0, -4, -9, -14, -23])

Vi_acc2 = np.array([180, 190, 200, 210, 220])
Fs_acc2 = np.array([0, -5, -9, -15, -24])

rho0 = 1.225
sigma = 0.57
S = 39.43
W_kgf = 10800
W1 = 11150 * 9.81
W2 = 10450 * 9.81

kt_to_ms = 0.514444
V_tas1 = (Vi_stab1 / np.sqrt(sigma)) * kt_to_ms
V_tas2 = (Vi_stab2 / np.sqrt(sigma)) * kt_to_ms

q1 = 0.5 * (rho0 * sigma) * V_tas1**2
q2 = 0.5 * (rho0 * sigma) * V_tas2**2

CL1 = W1 / (q1 * S)
CL2 = W2 / (q2 * S)

# A) Manche Fixo (Gradiente d(delta_e)/dCL)
slope_fixo1, intercept_fixo1 = np.polyfit(CL1, de_stab1, 1)
slope_fixo2, intercept_fixo2 = np.polyfit(CL2, de_stab2, 1)

hn_fixo = (h_cg1 * slope_fixo1 - h_cg2 * slope_fixo2) / (slope_fixo1 - slope_fixo2)
SM_fixo1 = hn_fixo - h_cg1
SM_fixo2 = hn_fixo - h_cg2

# B) Manche Livre (Gradiente dFs/dV)
slope_livre1, intercept_livre1 = np.polyfit(Vi_acc1, Fs_acc1, 1)
slope_livre2, intercept_livre2 = np.polyfit(Vi_acc2, Fs_acc2, 1)

hn_livre = (slope_livre1 * h_cg2 - slope_livre2 * h_cg1) / (slope_livre1 - slope_livre2)
SM_livre1 = hn_livre - h_cg1
SM_livre2 = hn_livre - h_cg2

print("-" * 40)
print("RESULTADOS DE ESTABILIDADE ESTÁTICA")
print("-" * 40)
print(f"Ponto Neutro Manche FIXO (hn):  {hn_fixo*100:.2f} % CMA")
print(f"Margem Estática Fixo (Voo 1):   {SM_fixo1*100:.2f} %")
print(f"Margem Estática Fixo (Voo 2):   {SM_fixo2*100:.2f} %")
print("-" * 40)
print(f"Ponto Neutro Manche LIVRE (hn'): {hn_livre*100:.2f} % CMA")
print(f"Margem Estática Livre (Voo 1):   {SM_livre1*100:.2f} %")
print(f"Margem Estática Livre (Voo 2):   {SM_livre2*100:.2f} %")
print("-" * 40)
print(f"Gradiente Manche Fixo Voo 1: {slope_fixo1:.2f} cm/CL")
print(f"Gradiente Manche Fixo Voo 2: {slope_fixo2:.2f} cm/CL")
print(f"Gradiente Força Voo 1: {slope_livre1:.3f} kgf/kt")
print(f"Gradiente Força Voo 2: {slope_livre2:.3f} kgf/kt")
print("-" * 40)

palette = sns.color_palette("Set2", n_colors=2)
plt.rcParams['font.family'] = 'Segoe UI'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(CL1, de_stab1, 'o-', label=f'Voo 1 (CG {h_cg1*100:.2f}%)', color=palette[0])
ax1.plot(CL2, de_stab2, 's-', label=f'Voo 2 (CG {h_cg2*100:.2f}%)', color=palette[1])
x_fit = np.linspace(0.4, 1.1, 10)
# ax1.plot(x_fit, slope_fixo1*x_fit + intercept_fixo1, '--', alpha=0.5, color=palette[0])
# ax1.plot(x_fit, slope_fixo2*x_fit + intercept_fixo2, '--', alpha=0.5, color=palette[1])

ax1.set_title('Estabilidade Estática Longitudinal (Manche Fixo)')
ax1.set_xlabel('Coeficiente de Sustentação ($C_L$)', fontsize=13)
ax1.set_ylabel('Posição do Manche $\delta_e$ (cm)', fontsize=13)
ax1.invert_xaxis()
ax1.invert_yaxis()
ax1.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax1.grid(axis='x')
ax1.legend()
sns.despine(ax=ax1)
for spine in ax1.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

ax2.plot(Vi_acc1, Fs_acc1, 'o-', label=f'Voo 1 (CG {h_cg1*100:.2f}%)', color=palette[0])
ax2.plot(Vi_acc2, Fs_acc2, 's-', label=f'Voo 2 (CG {h_cg2*100:.2f}%)', color=palette[1])
v_fit = np.linspace(180, 220, 10)
# ax2.plot(v_fit, slope_livre1*v_fit + intercept_livre1, '--', alpha=0.5, color=palette[0])
# ax2.plot(v_fit, slope_livre2*v_fit + intercept_livre2, '--', alpha=0.5, color=palette[1])

ax2.set_title('Estabilidade Estática Longitudinal (Manche Livre)')
ax2.set_xlabel('Velocidade Indicada ($V_{IAS}$) [kt]', fontsize=13)
ax2.set_ylabel('Força no Manche $F_s$ [kgf]', fontsize=13)
ax2.grid(axis='y', linestyle='--', linewidth=0.8, alpha=0.6)
ax2.grid(axis='x')
ax2.legend()
sns.despine(ax=ax2)
for spine in ax2.spines.values():
    spine.set_color('gray')
    spine.set_linewidth(1)

plt.tight_layout()
plt.savefig('estabilidade_longitudinal.png')
plt.show()