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

# Inicializar la simulación de Rebound
sim = rebound.Simulation()
sim.G = 4 * np.pi**2  # Unidades: AU^3 / (M_sun * yr^2)

# Añadir el Sol
sim.add(m=1.0)  # Masa del Sol en unidades solares

# Añadir Júpiter con masa total y elementos orbitales
sim.add(
    m=M_jup.value / M_sun.value,  # Masa de Júpiter en unidades solares
    a=0.01,  # Semieje mayor en AU
    e=0.4,  # Excentricidad
    f=np.pi  # Anomalía verdadera para comenzar en el apoastro
)

# Mover al centro de masa
sim.move_to_com()

# Configurar el integrador
sim.integrator = "ias15"  # Integrador preciso
sim.dt = 0.001  # Paso de integración pequeño

# Integrar la simulación para graficar la evolución de Júpiter
tiempo_total = 2.0  # Integrar durante 2 años
N_steps = 1000  # Número de pasos de la simulación

x_trayectoria = []  # Almacenar la trayectoria de Júpiter en x
y_trayectoria = []  # Almacenar la trayectoria de Júpiter en y

for _ in range(N_steps):
    sim.integrate(sim.t + tiempo_total / N_steps)
    p = sim.particles[1]  # Júpiter es el segundo cuerpo
    x_trayectoria.append(p.x)
    y_trayectoria.append(p.y)

# Configuración de la figura para graficar la evolución
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_xlabel('Distancia en UA')
ax.set_ylabel('Distancia en UA')
ax.set_xlim(-0.02, 0.02)
ax.set_ylim(-0.02, 0.02)
ax.set_title('Evolución de la órbita de Júpiter')

# Graficar el Sol en el foco
ax.plot(0, 0, 'o', color='yellow', markersize=20, label='Sol')

# Graficar el límite de Roche
theta_roche = np.linspace(0, 2 * np.pi, 500)
x_roche = d_roche * np.cos(theta_roche)
y_roche = d_roche * np.sin(theta_roche)
ax.plot(x_roche, y_roche, '--', color='red', label='Límite de Roche')


# Graficar la posición inicial de Júpiter
ax.plot(x_trayectoria[0], y_trayectoria[0], 'o', color='orange', markersize=10, label='Júpiter (posición inicial)')

# Graficar la trayectoria de Júpiter
ax.plot(x_trayectoria, y_trayectoria, '-', color='blue', alpha=0.7, linewidth=1, label='Júpiter')

# Mostrar el gráfico
plt.legend()
plt.show()

orbit_jupiter = rebound.Orbit(
    sim.particles[1].x,  # Posición X de Júpiter
    sim.particles[1].y,  # Posición Y de Júpiter
    sim.particles[1].z,  # Posición Z de Júpiter
    sim.particles[1].vx, # Velocidad X de Júpiter
    sim.particles[1].vy, # Velocidad Y de Júpiter
    sim.particles[1].vz, # Velocidad Z de Júpiter
    primary_mass=sim.particles[0].m  # Masa del Sol en unidades de masa de Rebound
)

# Mostrar los elementos orbitales de Júpiter
print("\nElementos Orbitales de Júpiter:")
print(f"Semieje Mayor (a): {orbit_jupiter.a:.2e} AU")
print(f"Excentricidad (e): {orbit_jupiter.e:.2f}")
print(f"Inclinación (i): {orbit_jupiter.inc:.2f} grados")
print(f"Longitud del Nodo Ascendente (Omega): {orbit_jupiter.Omega:.2f} grados")
print(f"Argumento del Periapsis (omega): {orbit_jupiter.omega:.2f} grados")
print(f"Anomalía Verdadera (f): {orbit_jupiter.f:.2f} grados")



