import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# importing display from jupyther
from IPython.display import display
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm

df = pd.read_csv(
    'https://docs.google.com/spreadsheets/d/e/2PACX-1vSr7tU5tK8cvvR32yypE1PArgXhmNbFJ9bw8w6Sm2zQhyaMs27csoo-77vgFedHw9z25Ez3Qm-geKyU/pub?gid=1806272337&single=true&output=csv'
)

print(df.head(10))


W0 = 48280
T0 = 125

plt.figure(figsize=(12, 6))
plt.rcParams['font.family'] = 'Segoe UI'

sns.scatterplot(data=df, x = 'Passageiros (Máx)', y = 'MTOW (kg)', color = 'skyblue', s=80, edgecolor='gray', alpha=0.6)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("Número de passageiros", fontsize=20)
plt.ylabel("MTOW (kg)", fontsize=20)

# Formatação dos eixos (separador de milhar com ponto)
plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))

# Colocando os nomes dos aviões como rótulos

# escolhendo as posições dos rótulos
#['Embraer E-175', 'Mitsubishi SpaceJet M90', 'Sukhoi Superjet 100', 'COMAC ARJ21-700', 'ATR 72-600', 'Boeing 717-200', 'Airbus A220-100', 'Embraer E190-E2', 'DHC Dash 8-400', 'Bombardier CRJ-900']
# dicionário com os desvios em y para cada rótulo
desvios = {
    'Embraer E-175': [0, 15],
    'Mitsubishi SpaceJet M90': [0, 15],
    'Sukhoi Superjet 100': [0, -15],
    'COMAC ARJ21-700': [10, 15],
    'ATR 72-600': [0, 12],
    'Boeing 717-200': [0, -20],
    'Airbus A220-100': [0, -20],
    'Embraer E190-E2': [0, -15],
    'DHC Dash 8-400': [0, 12],
    'Bombardier CRJ-900': [0, 15]}

# pegando os desvios em x e y
desvios_em_x = [val[0] for val in desvios.values()]
desvios_em_y = [val[1] for val in desvios.values()]

for i in range(len(df)):
    plt.annotate(df['Nome curto'][i],
             (df['Passageiros (Máx)'][i], df['MTOW (kg)'][i]),
             xytext=(desvios_em_x[i], desvios_em_y[i]),
             textcoords='offset points',
             fontsize=14,
             ha='center', va='center')

# plotando a linha de tendência com numpy
import numpy as np
Z = np.polyfit(df['Passageiros (Máx)'], df['MTOW (kg)'], 1)
p = np.poly1d(Z)
plt.plot(df['Passageiros (Máx)'], p(df['Passageiros (Máx)']), color='red', label=f'Tendência: {Z[0]:.2f}x + {Z[1]:.2f}', linestyle=':', linewidth=0.5)

#plotando o ponto de projeto
plt.scatter(92, W0, color='red', s=80, edgecolor='red', alpha=0.6)
plt.annotate('PDM-26',
             (92, W0), xytext=(0, 15), textcoords='offset points',
             fontsize=14, ha='center', va='center') 

# Grade leve
plt.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
plt.xlim(75, 120)
# Estilo limpo
sns.despine()
plt.tight_layout()
ax = plt.gca()


# Cor e espessura dos eixos
ax.spines['bottom'].set_color('gray')
ax.spines['left'].set_color('gray')
ax.spines['top'].set_color('gray')
ax.spines['right'].set_color('gray')

ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
plt.savefig('mtow_vs_passengers.png', dpi=300, bbox_inches='tight')
plt.show()



###

coluna1 = 'Takeoff_distance(m)'
coluna2 = 'Landing_distance(m)'
df[coluna1] = pd.to_numeric(df[coluna1].astype(str).str.replace(',', '').str.strip(), errors='coerce')
df[coluna2] = pd.to_numeric(df[coluna2].astype(str).str.replace(',', '').str.strip(), errors='coerce')

plt.figure(figsize=(12, 6))
plt.rcParams['font.family'] = 'Segoe UI'


sns.scatterplot(data=df, y = coluna2, x = coluna1, color = 'skyblue', s=80, edgecolor='gray', alpha=0.6)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(coluna1.replace("_", " "), fontsize=20)
plt.ylabel(coluna2.replace("_", " "), fontsize=20)

# Formatação dos eixos (separador de milhar com ponto)
plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))

# escolhendo as posições dos rótulos
desvios = {
    'Embraer E-175': [0, 15],
    'Sukhoi Superjet 100': [0, -15],
    'COMAC ARJ21-700': [0, 15],
    'ATR 72-600': [5, 12],
    'Boeing 717-200': [0, -15],
    'Airbus A220-100': [0, -20],
    'Embraer E190-E2': [0, -15],
    'DHC Dash 8-400': [0, 12],
    'Bombardier CRJ-900': [15, 15]}


# pegando os desvios em x e y
desvios_em_x = [val[0] for val in desvios.values()]
desvios_em_y = [val[1] for val in desvios.values()]



for i in range(len(df)):
    plt.annotate(df['Nome curto'][i],
             (df[coluna1][i], df[coluna2][i]),
             xytext=(desvios_em_x[i], desvios_em_y[i]),
             textcoords='offset points',
             fontsize=14,
             ha='center', va='center')
    
TO_dist = 1720  # Distância de decolagem do PDM-26
LD_dist = 1110  # Distância de pouso do PDM-26

#plotando o ponto de projeto
plt.scatter(TO_dist, LD_dist, color='red', s=80, edgecolor='red', alpha=0.6)
plt.annotate('PDM-26',
             (TO_dist, LD_dist), xytext=(0, 15), textcoords='offset points',
             fontsize=14, ha='center', va='center')

# Grade leve
plt.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
# Estilo limpo
sns.despine()
plt.tight_layout()
ax = plt.gca()
plt.ylim(1000,)
plt.xlim(1600, 2000)
# Cor e espessura dos eixos
ax.spines['bottom'].set_color('gray')
ax.spines['left'].set_color('gray')
ax.spines['top'].set_color('gray')
ax.spines['right'].set_color('gray')

ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
plt.savefig('takeoff_vs_landing_distance_concorr.png', dpi=300, bbox_inches='tight')
plt.show()

############ empuxo com mtow

plt.figure(figsize=(12, 6))

plt.rcParams['font.family'] = 'Segoe UI'

sns.scatterplot(data=df, y = 'Empuxo (kN)', x = 'MTOW (kg)', color = 'skyblue', s=80, edgecolor='gray', alpha=0.6)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("MTOW (kg)", fontsize=20)
plt.ylabel("Empuxo Máximo (kN)", fontsize=20)

# Formatação dos eixos (separador de milhar com ponto)
plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}".replace(",", ".")))

# escolhendo as posições dos rótulos
desvios = {
    'Embraer E-175': [-15, 15],
    'Mitsubishi SpaceJet M90': [0, 15],
    'Sukhoi Superjet 100': [0, -15],
    'COMAC ARJ21-700': [0, 15],
    'ATR 72-600': [5, 12],
    'Boeing 717-200': [0, -15],
    'Airbus A220-100': [0, -20],
    'Embraer E190-E2': [0, -15],
    'DHC Dash 8-400': [0, 12],
    'Bombardier CRJ-900': [15, 15]}

# pegando os desvios em x e y
desvios_em_x = [val[0] for val in desvios.values()]
desvios_em_y = [val[1] for val in desvios.values()]



for i in range(len(df)):
    plt.annotate(df['Nome curto'][i],
             (df['MTOW (kg)'][i], df['Empuxo (kN)'][i]),
             xytext=(desvios_em_x[i], desvios_em_y[i]),
             textcoords='offset points',
             fontsize=14,
             ha='center', va='center')


#plotando o ponto de projeto
plt.scatter(W0, T0, color='red', s=80, edgecolor='red', alpha=0.6)
plt.annotate('PDM-26',
             (W0, T0), xytext=(0, 15), textcoords='offset points',
             fontsize=14, ha='center', va='center') 

# Grade leve
plt.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.6)
# Estilo limpo
sns.despine()
plt.tight_layout()
ax = plt.gca()

plt.ylim(120, 220)

# Cor e espessura dos eixos
ax.spines['bottom'].set_color('gray')
ax.spines['left'].set_color('gray')
ax.spines['top'].set_color('gray')
ax.spines['right'].set_color('gray')

ax.spines['bottom'].set_linewidth(1)
ax.spines['left'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
plt.savefig('empuxo_vs_mtow.png', dpi=300, bbox_inches='tight')
plt.show()



