import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import G, M_sun, M_jup, R_jup

# Función para calcular el límite de Roche
def roche_limit(M_star, M_planet, R_planet):
    return R_planet * (2 * M_star / M_planet)**(1/3)

# Cálculo del límite de Roche con unidades
d_roche = roche_limit(M_sun, M_jup, R_jup).to(u.au)

# Solicitar al usuario la distancia de afelio y la excentricidad
#Apoastro/Afelio: Punto de la órbita donde el planeta está más alejado del Sol.
r_apoastro    = 0.01* u.au # float(input("Ingresa la distancia de afelio de Júpiter al Sol en UA: ")) * u.au
excentricidad = 0.4 #float(input("Ingresa la excentricidad de la órbita de Júpiter (entre 0 y 1): "))

# Calcular el semieje mayor y menor
a = r_apoastro / (1 + excentricidad)  # Semieje mayor
b = a * np.sqrt(1 - excentricidad**2)  # Semieje menor

# Calcular la distancia del perihelio
#Perihelio: Punto de la órbita donde el planeta está más cercano al Sol.
r_periastro = a * (1 - excentricidad)

# Crear el ángulo theta para graficar la órbita
theta = np.linspace(0, 2 * np.pi, 500)

# Ecuaciones paramétricas de la elipse ajustadas para iniciar en el afelio
x_orbita = (a * np.cos(theta) + a * excentricidad).to(u.au).value  # Desplazamos el centro de la elipse
y_orbita = (b * np.sin(theta)).to(u.au).value

# Configuración de la figura
fig, ax = plt.subplots(figsize=(8,8))
ax.set_aspect('equal')
ax.set_xlabel('Distancia en UA')
ax.set_ylabel('Distancia en UA')
ax.set_title('Órbita de Júpiter y Límite de Roche')

# Graficar el Sol en el foco
ax.plot(0, 0, 'o', color='yellow', markersize=20, label='Sol')

# Graficar la órbita de Júpiter
ax.plot(x_orbita, y_orbita, color='blue', label='Órbita de Júpiter')

# Posición inicial de Júpiter en el afelio
ax.plot(r_apoastro.to(u.au).value, 0, 'o', color='orange', markersize=10, label='Júpiter en afelio')

# Graficar el límite de Roche
theta_roche = np.linspace(0, 2 * np.pi, 100)
x_roche = d_roche.value * np.cos(theta_roche)
y_roche = d_roche.value * np.sin(theta_roche)
ax.plot(x_roche, y_roche, '--', color='red', label='Límite de Roche')

# Ajustar los límites del gráfico
max_distance = r_apoastro.to(u.au).value * 1.1
ax.set_xlim(-max_distance, max_distance)
ax.set_ylim(-max_distance, max_distance)

# Agregar leyenda
ax.legend()

# Calcular el período orbital utilizando la tercera ley de Kepler
T = 2 * np.pi * np.sqrt(a**3 / (G * M_sun))
T = T.to(u.year)

# Mostrar información
print(f"\nEl período orbital del planeta es {T.value:.2e} años.")
if T.value < 1:
    T_dias = T.to(u.day)
    print(f"Que equivale a {T_dias.value:.2f} días.")
print(f"Apoastro: {r_apoastro.to(u.au).value:.2e} UA.")
print(f"Periastro: {r_periastro.to(u.au).value:.2e} UA.")
print(f"Límite de Roche: {d_roche.value:.2e} UA.")
print(f"Excentricidad: {excentricidad}")

# Mostrar el gráfico
plt.show()
