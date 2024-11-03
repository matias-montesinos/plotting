import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from scipy.signal import savgol_filter  # Importar savgol_filter
from scipy.interpolate import UnivariateSpline


# Rutas donde están almacenados los datos
path_1MJ_YESfeel_lower = "/home/matias/Simulations/mi_fargo3d/outputs/flyby2d_1MJ_YESfeel_lower/"
path = path_1MJ_YESfeel_lower

# Leer variables.par
variables_par = np.genfromtxt(path + "/variables.par", dtype={'names': ("parametros", "valores"), 'formats': ("|S30", "|S300")}).tolist()
parametros_par, valores_par = [], []
for posicion in variables_par:
    parametros_par.append(posicion[0].decode("utf-8"))
    valores_par.append(posicion[1].decode("utf-8"))

def P(parametro):
    return valores_par[parametros_par.index(parametro)] 

# Parámetros de la simulación
Ninter = int(P("NINTERM"))
DT = float(P("DT"))
NX = int(P("NX"))
NY = int(P("NY"))

# Obtener Rmax
try:
    Rmax = float(P("RMAX"))
except:
    domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]
    Rmax = domain_y[-1]

def obtener_momento_angular_y_distancia(path):
    """
    Calcula el momento angular específico y la distancia radial del planeta.
    """
    planet_file = path + "planet0.dat"
    planet_data = np.genfromtxt(planet_file)
    
    x_coords = planet_data[:, 1]
    y_coords = planet_data[:, 2]
    v_azimutal = planet_data[:, 4]
    
    r_distances = np.sqrt(x_coords**2 + y_coords**2)
    momento_angular_especifico_planeta = r_distances * v_azimutal
    
    return momento_angular_especifico_planeta, r_distances

# Obtener los datos
momento_angular, r_distances = obtener_momento_angular_y_distancia(path)

# Crear array de snapshots
snapshots = np.arange(len(momento_angular))

# Identificar la condición donde el planeta está dentro del disco
inside_disk = r_distances < Rmax

# Crear segmentos para LineCollection
points = np.array([snapshots, momento_angular]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Crear un array de colores basado en la condición
colors = np.where(inside_disk[:-1], 'red', 'blue')

# Crear LineCollection
lc = LineCollection(segments, colors=colors, linewidths=2)

# Crear la figura y los ejes
fig, ax = plt.subplots()

# Añadir la LineCollection a los ejes
ax.add_collection(lc)

# Ajustar los límites de los ejes
ax.set_xlim(snapshots.min(), snapshots.max())
ax.set_ylim(momento_angular.min(), momento_angular.max())

# Añadir etiquetas y título
ax.set_xlabel("Número de Snapshot")
ax.set_ylabel("Momento Angular Específico [UA² / año]")
ax.set_title("Evolución del Momento Angular Específico del Planeta")

# Añadir leyenda personalizada
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], color='red', lw=2, label='Dentro del disco'),
                   Line2D([0], [0], color='blue', lw=2, label='Fuera del disco')]
ax.legend(handles=legend_elements)
ax.grid(True)
plt.show()

# Paso 1: Calcular el array de tiempos para cada snapshot
delta_t = Ninter * DT  # Tiempo entre snapshots
time_array = snapshots #* delta_t  # Tiempo en cada snapshot

# Paso 2: Calcular la tasa de cambio del momento angular (derivada temporal)
delta_momentum = np.diff(momento_angular)
delta_time = np.diff(time_array)
rate_of_change = delta_momentum / delta_time  # Tasa de cambio del momento angular

# Paso 3: Preparar los datos para el gráfico
# Como diff reduce la longitud del array en 1, ajustamos los arrays correspondientes
time_midpoints = time_array[:-1] + delta_time / 2  # Puntos medios en el tiempo para las derivadas
inside_disk_midpoints = inside_disk[:-1]  # Usamos el estado 'dentro/fuera' del snapshot inicial de cada intervalo

# Paso 4: Crear segmentos para LineCollection
points = np.array([time_midpoints, rate_of_change]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Paso 5: Definir colores basados en la condición
colors = np.where(inside_disk_midpoints[:-1], 'red', 'blue')

# Paso 6: Crear LineCollection para el gráfico
lc = LineCollection(segments, colors=colors, linewidths=2)


# Paso 7: Graficar
# Define la variable de control al inicio
use_log_scale = False  # Cambia a False si no deseas escala logarítmica

fig, ax = plt.subplots()

ax.add_collection(lc)
ax.set_ylim(rate_of_change.min(), rate_of_change.max())

# Configurar el eje x según la variable use_log_scale
if use_log_scale:
    # Asegurarse de que no haya valores de tiempo cero o negativos
    valid_times = time_midpoints > 0
    ax.set_xscale('log')  # Establece escala logarítmica en el eje x
    ax.set_xlim(time_midpoints[valid_times].min(), time_midpoints.max())
else:
    ax.set_xlim(time_midpoints.min(), time_midpoints.max())

ax.set_xlabel("Tiempo [unidades de tiempo]")
ax.set_ylabel("Tasa de cambio del Momento Angular [UA² / año²]")
ax.set_title("Tasa de Cambio del Momento Angular del Planeta")

# Añadir leyenda personalizada
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='Dentro del disco'),
    Line2D([0], [0], color='blue', lw=2, label='Fuera del disco')
]
ax.legend(handles=legend_elements)

ax.grid(True)
plt.show()


# Parámetros del filtro de Savitzky-Golay
window_length = 11  # Debe ser un número entero positivo impar
polyorder = 3       # Orden del polinomio ajustado

# Ajustar window_length si es necesario
if window_length >= len(momento_angular):
    window_length = len(momento_angular) - 1 if len(momento_angular) % 2 == 1 else len(momento_angular) - 2

# Aplicar el filtro para suavizar los datos y calcular la derivada
delta_t = Ninter * DT  # Ya definido en tu código

# Suavizar el momento angular
momento_angular_suave = savgol_filter(momento_angular, window_length, polyorder)

# Calcular la derivada suavizada
derivada_momento_angular = savgol_filter(momento_angular, window_length, polyorder, deriv=1, delta=delta_t)

# Ajustar spline
spline = UnivariateSpline(time_array, momento_angular, s=0, k=3)
# Calcular derivada
#derivada_momento_angular = spline.derivative()(time_array)


# Crear array de tiempos
time_array = snapshots * delta_t

# Asegurarse de que inside_disk tenga la misma longitud que time_array
inside_disk = r_distances < Rmax
# Puntos para los segmentos
points = np.array([time_array, derivada_momento_angular]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Colores basados en la condición
colors = np.where(inside_disk[:-1], 'red', 'blue')
# Crear LineCollection
lc = LineCollection(segments, colors=colors, linewidths=2)

# Graficar
fig, ax = plt.subplots()

ax.add_collection(lc)
ax.set_xlim(time_array.min(), time_array.max())
ax.set_ylim(derivada_momento_angular.min(), derivada_momento_angular.max())

ax.set_xlabel("Tiempo [unidades de tiempo]")
ax.set_ylabel("Tasa de cambio del Momento Angular [UA² / año²]")
ax.set_title("Tasa de Cambio Suavizada del Momento Angular del Planeta")

# Añadir leyenda personalizada
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='Dentro del disco'),
    Line2D([0], [0], color='blue', lw=2, label='Fuera del disco')
]
ax.legend(handles=legend_elements)

ax.grid(True)
plt.show()



###
# Ajustar window_length si es necesario
if window_length >= len(momento_angular):
    window_length = len(momento_angular) - 1 if len(momento_angular) % 2 == 1 else len(momento_angular) - 2
    if window_length < 3:
        window_length = 3  # Valor mínimo aceptable

# Aplicar el filtro para suavizar el momento angular
momento_angular_suave = savgol_filter(momento_angular, window_length, polyorder)

# Paso 3: Crear segmentos para LineCollection
# Puntos para los segmentos
points = np.array([time_array, momento_angular_suave]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Colores basados en la condición
colors = np.where(inside_disk[:-1], 'red', 'blue')

# Paso 4: Graficar el momento angular suavizado con colores dentro/fuera del disco
fig, ax = plt.subplots()

# Crear LineCollection
lc = LineCollection(segments, colors=colors, linewidths=2)

# Añadir la colección al gráfico
ax.add_collection(lc)

# Ajustar los límites de los ejes
ax.set_xlim(time_array.min(), time_array.max())
ax.set_ylim(momento_angular_suave.min(), momento_angular_suave.max())

# Añadir etiquetas y título
ax.set_xlabel("Tiempo [unidades de tiempo]")
ax.set_ylabel("Momento Angular Específico [UA² / año]")
ax.set_title("Evolución Suavizada del Momento Angular Específico del Planeta")

# Añadir leyenda personalizada
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='Dentro del disco'),
    Line2D([0], [0], color='blue', lw=2, label='Fuera del disco')
]
ax.legend(handles=legend_elements)

# Mostrar la cuadrícula
ax.grid(True)

# Mostrar el gráfico
plt.show()




