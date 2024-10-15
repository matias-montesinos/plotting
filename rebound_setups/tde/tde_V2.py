
import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
from astropy.constants import G, M_sun, M_jup, R_jup
import os

# Función para calcular el límite de Roche
def roche_limit(M_star, M_planet, R_planet):
    return R_planet * (2 * M_star / M_planet)**(1/3)

# Cálculo del límite de Roche con unidades
d_roche = roche_limit(M_sun, M_jup, R_jup).to(u.au)

# Definir los parámetros de la órbita
r_apoastro = 0.01 * u.au         # Distancia de apoastro en AU
excentricidad = 0.4               # Excentricidad de la órbita

# Calcular el semieje mayor
a = r_apoastro / (1 + excentricidad)  # Semieje mayor

# Calcular el período orbital utilizando la tercera ley de Kepler
T = 2 * np.pi * np.sqrt(a**3 / (G * M_sun))
T = T.to(u.year)

# Mostrar información preliminar
print(f"\nSemieje mayor (a): {a.to(u.au).value:.2e} UA")
print(f"Límite de Roche: {d_roche.value:.2e} UA")
print(f"Excentricidad: {excentricidad}")
print(f"Período orbital: {T.value:.2e} años.")

if T.value < 1:
    T_dias = T.to(u.day)
    print(f"Que equivale a {T_dias.value:.2f} días.")

# Inicializar la simulación de Rebound
sim = rebound.Simulation()
sim.G = 4 * np.pi**2  # Unidades: AU^3 / (M_sun * yr^2)

# Añadir el Sol
sim.add(m=1.0)  # Masa del Sol en unidades solares

# Añadir Júpiter en el apoastro con anomalía verdadera f=pi
sim.add(
    m=M_jup.value / M_sun.value,  # Masa de Júpiter en unidades solares
    a=a.to(u.au).value,            # Semieje mayor en AU
    e=excentricidad,               # Excentricidad
    f=np.pi                         # Anomalía verdadera para apoastro
)

# Verificar posición inicial de Júpiter
jupiter = sim.particles[1]
print(f"\nPosición inicial de Júpiter: x = {jupiter.x:.2e} AU, y = {jupiter.y:.2e} AU")
print(f"Velocidad inicial de Júpiter: vx = {jupiter.vx:.2e} AU/yr, vy = {jupiter.vy:.2e} AU/yr")

# Mover al centro de masa
sim.move_to_com()

# Configuración de la figura
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlabel('Distancia en UA')
ax.set_ylabel('Distancia en UA')
ax.set_title('Órbita de Júpiter y Límite de Roche con Rebound')

# Graficar el Sol en el foco
ax.plot(0, 0, 'o', color='yellow', markersize=20, label='Star')

# Integrar la simulación para graficar la órbita de Júpiter
tiempo_total = T.value
N_steps = 1000
tiempos = np.linspace(0, tiempo_total, N_steps)
x_rebound = np.zeros(N_steps)
y_rebound = np.zeros(N_steps)

for i, t in enumerate(tiempos):
    sim.integrate(t)
    p = sim.particles[1]  # Júpiter es el segundo cuerpo
    x_rebound[i] = p.x
    y_rebound[i] = p.y

# Invertir el eje X para que el apoastro esté en +x
x_orbita = -x_rebound
y_orbita = y_rebound

# Graficar la órbita de Júpiter
ax.plot(x_orbita, y_orbita, color='blue', label='Planet orbit')

# Posición inicial de Júpiter en el apoastro
ax.plot(r_apoastro.to(u.au).value, 0, 'o', color='orange', markersize=10, label='Júpiter en apoastro')

# Graficar el límite de Roche
theta_roche = np.linspace(0, 2 * np.pi, 500)
x_roche = d_roche.value * np.cos(theta_roche)
y_roche = d_roche.value * np.sin(theta_roche)
ax.plot(x_roche, y_roche, '--', color='red', label='Roche limit')

# Ajustar los límites del gráfico
max_distance = r_apoastro.to(u.au).value * 1.1
ax.set_xlim(-max_distance, max_distance)
ax.set_ylim(-max_distance, max_distance)

# Agregar leyenda
ax.legend()

# Mostrar el gráfico
plt.show()

# Utilizar rebound.OrbitPlot para visualizar la órbita
#op = rebound.OrbitPlot(sim, color=True, unitlabel="[AU]")
#plt.show() 

# Mostrar información adicional
print(f"\nEl período orbital del planeta es {T.value:.2e} años.")
if T.value < 1:
    T_dias = T.to(u.day)
    print(f"Que equivale a {T_dias.value:.2f} días.")
print(f"Apoastro: {r_apoastro.to(u.au).value:.2e} UA.")
print(f"Límite de Roche: {d_roche.value:.2e} UA.")
print(f"Excentricidad: {excentricidad}")


from matplotlib.animation import FuncAnimation
# Configuración de la figura para la animación
# Configuración de la figura para la animación
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlabel('X [AU]')
ax.set_ylabel('Y [AU]')
ax.set_xlim(-0.02, 0.02)
ax.set_ylim(-0.02, 0.02)
ax.set_title('Orbit evolution')

# Graficar el Sol en el foco
ax.plot(0, 0, 'o', color='yellow', markersize=20, label='Star')

#grafico limite Roche
ax.plot(x_roche, y_roche, '--', color='red', label='Roche limit')


# Inicializar la línea de la órbita de Júpiter
ojupiter, = ax.plot([], [], 'o', color='orange', markersize=10, label='Planet')
linea_orbita, = ax.plot([], [], '-', color='blue', lw=1.5, label='Planet\'s orbit')

# Inicializar los datos de la trayectoria
xdata, ydata = [], []

# Función de inicialización para la animación
def init():
    xdata.clear()
    ydata.clear()
    # Posición inicial de Júpiter en el apoastro
    x_apoastro = r_apoastro.to(u.au).value
    y_apoastro = 0.0
    xdata.append(x_apoastro)
    ydata.append(y_apoastro)
    ojupiter.set_data([x_apoastro], [y_apoastro])  # Colocar a Júpiter en la posición inicial correcta
    linea_orbita.set_data([], [])  # Iniciar la línea vacía
    return ojupiter, linea_orbita

# Función de actualización para cada cuadro de la animación
def update(frame):
    x = -x_rebound[frame]  # Invertir el eje X para mantener la coherencia con el gráfico estático
    y = -y_rebound[frame]
    xdata.append(x)
    ydata.append(y)
    ojupiter.set_data([x], [y])  # Pasar valores como listas para que set_data() funcione correctamente
    linea_orbita.set_data(xdata, ydata)
    return ojupiter, linea_orbita


# Crear la animación
ani = FuncAnimation(fig, update, frames=range(0, N_steps, 10), interval=50 ,init_func=init, blit=True)

output_dir = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d/"

# Crear el directorio si no existe
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_file = 'planet_orbit_with_gas_animation.mp4'
ani.save(output_file, writer='ffmpeg', fps=10)

# Mostrar la animación
plt.legend()
plt.show()