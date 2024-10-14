import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
from astropy.constants import G, M_sun, M_jup, R_jup

# Función para calcular el límite de Roche
def roche_limit(M_star, M_planet, R_planet):
    return R_planet * (2 * M_star / M_planet)**(1/3)

# Cálculo del límite de Roche con unidades
d_roche = roche_limit(M_sun, M_jup, R_jup).to(u.au)

# Definir los parámetros de la órbita
r_apoastro = 0.01 * u.au         # Distancia de apoastro en AU
excentricidad = 0.4               # Excentricidad de la órbita

# Calcular el semieje mayor y menor
a = r_apoastro / (1 + excentricidad)  # Semieje mayor
b = a * np.sqrt(1 - excentricidad**2)  # Semieje menor

# Calcular la distancia del periastro
r_periastro = a * (1 - excentricidad)

# Calcular el período orbital utilizando la tercera ley de Kepler
T = 2 * np.pi * np.sqrt(a**3 / (G * M_sun))
T = T.to(u.year)

# Mostrar información preliminar
print(f"\nSemieje mayor (a): {a.to(u.au).value:.2e} UA")
print(f"Periastro: {r_periastro.to(u.au).value:.2e} UA")
print(f"Apoastro: {r_apoastro.to(u.au).value:.2e} UA")
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

# Configurar el integrador
sim.integrator = "ias15"  # Integrador preciso
sim.dt = 0.001 * T.value  # Paso de integración pequeño

# Definir el tiempo de integración
tiempo_total = T.value  # Integrar durante un período orbital

# Definir el número de puntos para la trayectoria
N_steps = 1000
tiempos = np.linspace(0, tiempo_total, N_steps)
x_rebound = np.zeros(N_steps)
y_rebound = np.zeros(N_steps)

# Integrar y registrar posiciones
for i, t in enumerate(tiempos):
    sim.integrate(t)
    p = sim.particles[1]  # Júpiter es el segundo cuerpo
    x_rebound[i] = p.x
    y_rebound[i] = p.y

# Invertir el eje X para que el apoastro esté en +x
x_orbita = -x_rebound
y_orbita = y_rebound

# Configuración de la figura
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlabel('Distancia en UA')
ax.set_ylabel('Distancia en UA')
ax.set_title('Órbita de Júpiter y Límite de Roche con Rebound')

# Graficar el Sol en el foco
ax.plot(0, 0, 'o', color='yellow', markersize=20, label='Sol')

# Graficar la órbita de Júpiter
ax.plot(x_orbita, y_orbita, color='blue', label='Órbita de Júpiter')

# Posición inicial de Júpiter en el apoastro
ax.plot(r_apoastro.to(u.au).value, 0, 'o', color='orange', markersize=10, label='Júpiter en apoastro')

# Graficar el límite de Roche
theta_roche = np.linspace(0, 2 * np.pi, 500)
x_roche = d_roche.value * np.cos(theta_roche)
y_roche = d_roche.value * np.sin(theta_roche)
ax.plot(x_roche, y_roche, '--', color='red', label='Límite de Roche')

# Ajustar los límites del gráfico
max_distance = r_apoastro.to(u.au).value * 1.1
ax.set_xlim(-max_distance, max_distance)
ax.set_ylim(-max_distance, max_distance)

# Agregar leyenda
ax.legend()

# Mostrar el gráfico
plt.show()

# Mostrar información adicional
print(f"\nEl período orbital del planeta es {T.value:.2e} años.")
if T.value < 1:
    T_dias = T.to(u.day)
    print(f"Que equivale a {T_dias.value:.2f} días.")
print(f"Apoastro: {r_apoastro.to(u.au).value:.2e} UA.")
print(f"Periastro: {r_periastro.to(u.au).value:.2e} UA.")
print(f"Límite de Roche: {d_roche.value:.2e} UA.")
print(f"Excentricidad: {excentricidad}")


# Calcular las órbitas
orbits = sim.orbits()

orbit_jupiter = rebound.Orbit(
    sim.particles[1].x,  # Posición X de Júpiter
    sim.particles[1].y,  # Posición Y de Júpiter
    sim.particles[1].z,  # Posición Z de Júpiter
    sim.particles[1].vx, # Velocidad X de Júpiter
    sim.particles[1].vy, # Velocidad Y de Júpiter
    sim.particles[1].vz, # Velocidad Z de Júpiter
    primary_mass=sim.particles[0].m  # Masa del Sol en unidades de masa de Rebound
)
print("\nElementos Orbitales de Júpiter:")
print(f"Semieje Mayor (a): {orbit_jupiter.a:.2e} AU")
print(f"Excentricidad (e): {orbit_jupiter.e:.2f}")
print(f"Inclinación (i): {orbit_jupiter.inc:.2f} grados")
print(f"Longitud del Nodo Ascendente (Omega): {orbit_jupiter.Omega:.2f} grados")
print(f"Argumento del Periapsis (omega): {orbit_jupiter.omega:.2f} grados")
print(f"Anomalía Verdadera (f): {orbit_jupiter.f:.2f} grados")

# Después de configurar y ejecutar la simulación
op = rebound.OrbitPlot(sim)
