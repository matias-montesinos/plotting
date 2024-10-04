import numpy as np
import matplotlib.pyplot as plt

# Ruta donde está almacenado el archivo planet0.dat
path = "/Users/matias/Simulations/mi_fargo3d/outputs/flyby2d_1MJ_YESfeel_lower/"
planet_file = path + "planet0.dat"

# Leer el archivo planet0.dat (ajusta el delimitador si es necesario)
planet_data = np.genfromtxt(planet_file)

# Extraer las coordenadas x e y
x_coords = planet_data[:, 1]  # Segunda columna: coordenada x
y_coords = planet_data[:, 2]  # Tercera columna: coordenada y

# Graficar todas las coordenadas
plt.figure(figsize=(8, 8))

# Graficar todos los puntos intermedios en negro
plt.plot(x_coords, y_coords, color='black', linestyle='-', marker='o', markersize=5, label='Órbita del planeta')

# Marcar la posición inicial en azul
plt.scatter(x_coords[0], y_coords[0], color='blue', s=100, label='Posición inicial')

# Marcar la posición final en rojo
plt.scatter(x_coords[-1], y_coords[-1], color='red', s=100, label='Posición final')

# Agregar el sol en el centro (0,0) en amarillo
plt.scatter(0, 0, color='yellow', s=200, edgecolor='orange', label='Star', zorder=5)

# Configurar etiquetas y leyenda
plt.xlabel('X [AU]', fontsize=14)
plt.ylabel('Y [AU]', fontsize=14)
plt.title('Planet orbit', fontsize=16)
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Mantener proporciones

# Mostrar el gráfico
plt.show()
