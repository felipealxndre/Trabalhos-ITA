import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

file_path = 'C:\\Users\\fealp\\OneDrive\\Documentos\\ITA\\Trabalhos\\MVO-32\\REV\\v5_20251008_085011_105m.csv'
df = pd.read_csv(file_path)

df['Timestamp'] = pd.to_datetime(df['Timestamp'], format='%d-%b-%Y %H:%M:%S.%f')
df['Time_s'] = (df['Timestamp'] - df['Timestamp'].iloc[0]).dt.total_seconds()


t_start = 0
t_end = 5000

mask = (df['Time_s'] >= t_start) & (df['Time_s'] <= t_end)
df_short = df.loc[mask].copy()

time = df_short['Time_s'].values
q = df_short['gyroY(rad/s)'].values
pitch = df_short['Pitch(rads)'].values

peaks, _ = find_peaks(q, prominence=0.02)
valleys, _ = find_peaks(-q, prominence=0.02)


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)

ax1.plot(time, q, 'b-', linewidth=1.5, label='Taxa de Arfagem ($q$)')
ax1.plot(time[peaks], q[peaks], 'rx', markersize=10, label='Picos')
ax1.plot(time[valleys], q[valleys], 'go', markersize=8, label='Vales')
ax1.set_ylabel('Taxa de Arfagem $q$ [rad/s]', fontsize=12)
ax1.axhline(0, color='black', linewidth=0.8, linestyle='--')
ax1.legend()
ax1.grid(True, which='both', linestyle='--', alpha=0.7)
ax1.set_title(f'Análise de Curto Período (Intervalo {t_start}-{t_end}s)', fontsize=14)
ax1.minorticks_on()

ax2.plot(time, pitch, 'm-', linewidth=1.2, label='Pitch (rad)')
ax2.set_xlabel('Tempo [s]', fontsize=12)
ax2.set_ylabel('Pitch [rad]', fontsize=12)
ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
ax2.legend()
ax2.grid(True, which='both', linestyle='--', alpha=0.7)
ax2.minorticks_on()

plt.tight_layout()
plt.show()