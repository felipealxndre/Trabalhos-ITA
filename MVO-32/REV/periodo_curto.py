import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# 1. CARREGAMENTO E PREPARAÇÃO
file_path = 'C:\\Users\\fealp\\OneDrive\\Documentos\\ITA\\Trabalhos\\MVO-32\\REV\\v5_20251008_085011_105m.csv'
df = pd.read_csv(file_path)

# Ajustar Tempo
df['Timestamp'] = pd.to_datetime(df['Timestamp'], format='%d-%b-%Y %H:%M:%S.%f')
df['Time_s'] = (df['Timestamp'] - df['Timestamp'].iloc[0]).dt.total_seconds()

# 2. RECORTE DO TEMPO (Zoom na manobra)
# Intervalo identificado anteriormente como provável Curto Período
t_start = 0
t_end = 5000

mask = (df['Time_s'] >= t_start) & (df['Time_s'] <= t_end)
df_short = df.loc[mask].copy()

# Variável de Interesse: Taxa de Arfagem q (Giroscópio X)
time = df_short['Time_s'].values
q = df_short['gyroY(rad/s)'].values
pitch = df_short['Pitch(rads)'].values

# 3. ANÁLISE DE PICOS (Identificação da Dinâmica)
# Encontrar picos positivos e negativos (vales) para medir o período
# Ajuste 'prominence' se necessário para filtrar ruído
peaks, _ = find_peaks(q, prominence=0.02)
valleys, _ = find_peaks(-q, prominence=0.02)

# 4. PLOTAGEM

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

# Subplot 1: Taxa de arfagem q
ax1.plot(time, q, 'b-', linewidth=1.5, label='Taxa de Arfagem ($q$)')
ax1.plot(time[peaks], q[peaks], 'rx', markersize=10, label='Picos')
ax1.plot(time[valleys], q[valleys], 'go', markersize=8, label='Vales')
ax1.set_ylabel('Taxa de Arfagem $q$ [rad/s]', fontsize=12)
ax1.axhline(0, color='black', linewidth=0.8, linestyle='--')
ax1.legend()
ax1.grid(True, which='both', linestyle='--', alpha=0.7)
ax1.set_title(f'Análise de Curto Período (Intervalo {t_start}-{t_end}s)', fontsize=14)
ax1.minorticks_on()

# Subplot 2: Pitch
ax2.plot(time, pitch, 'm-', linewidth=1.2, label='Pitch (rad)')
ax2.set_xlabel('Tempo [s]', fontsize=12)
ax2.set_ylabel('Pitch [rad]', fontsize=12)
ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
ax2.legend()
ax2.grid(True, which='both', linestyle='--', alpha=0.7)
ax2.minorticks_on()

plt.tight_layout()
plt.show()

# 5. CÁLCULO DOS PARÂMETROS
print("--- RESULTADOS DA ANÁLISE ---")
if len(peaks) >= 2:
    # Calcular período médio entre picos consecutivos
    diff_peaks = np.diff(time[peaks])
    Td_medio = np.mean(diff_peaks)
    wn = (2 * np.pi) / Td_medio
    
    print(f"Picos detectados nos instantes: {time[peaks]}")
    print(f"Período Amortecido Médio (Td): {Td_medio:.3f} s")
    print(f"Frequência Natural Estimada (wn): {wn:.3f} rad/s")
    
    # Estimativa de Amortecimento (Decremento Logarítmico)
    # Usando o primeiro e segundo pico (se houver decaimento claro)
    if q[peaks][0] > q[peaks][1]:
        delta = np.log(q[peaks][0] / q[peaks][1])
        zeta = delta / np.sqrt(4 * np.pi**2 + delta**2)
        print(f"Amortecimento Estimado (Decremento Log): {zeta:.3f}")
    else:
        print("Amortecimento: Não foi possível calcular (picos não decrescentes ou ruído).")
        
else:
    print("Não foram detectados picos suficientes para cálculo automático.")
    print("Tente ajustar o intervalo de tempo ou a proeminência dos picos.")

