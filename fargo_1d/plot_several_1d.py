import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re

path = "/home/matias/Simulations/fargo3d/outputs/accretion/"

# Leer el archivo de parámetros
variables_par = np.genfromtxt(path + "/variables.par", dtype={'names': ("parametros", "valores"), 'formats': ("|S30", "|S300")}).tolist()
parametros_par, valores_par = [], []
for posicion in variables_par:
    parametros_par.append(posicion[0].decode("utf-8"))
    valores_par.append(posicion[1].decode("utf-8"))

def P(parametro):
    return valores_par[parametros_par.index(parametro)]

# Identificar el número de outputs disponibles
pattern = re.compile(r'gasdens(\d+)\.dat')
files = glob.glob(path + "gasdens*.dat")
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
Ntot = len(valid_files)
print(f"Número total de outputs disponibles: {Ntot}")

# Preguntar al usuario el intervalo entre outputs a plotear
intervalo = int(input(f"Ingresa el intervalo entre outputs a plotear (ej. 2 para cada 2): "))

# Crear lista de outputs a plotear
outputs_a_plotear = list(range(0, Ntot, intervalo))
print(f"Outputs a plotear: {outputs_a_plotear}")

# Datos de la grilla
Rmin = float(P("YMIN"))
Rmax = float(P("YMAX"))

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]  # Cerrar el corchete

NX = int(P("NX")) 
NY = int(P("NY")) 
dr = (Rmax - Rmin) / NY


rmin = domain_y
rmed = 0.5 * (rmin[1:] + rmin[:-1])

# Coordenada r
r_values = rmed
print(r_values.shape)

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out

# Función para cargar vel radial
def cargar_vrad(output_number, path):
    v_out = np.fromfile(path + f"gasvy{output_number}.dat").reshape(NY, NX)
    return v_out

# Función para cargar vel azimuthal
def cargar_vazi(output_number, path):
    v_out = np.fromfile(path + f"gasvx{output_number}.dat").reshape(NY, NX)
    return v_out

# Configuración de los gráficos
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

# Plotear todos los outputs seleccionados en el mismo gráfico
for snapshot in outputs_a_plotear:
    # Cargar los datos de cada output
    dens_out = cargar_densidad(snapshot, path)
    vrad_out = cargar_vrad(snapshot, path)
    vazi_out = cargar_vazi(snapshot, path)
    #print(vrad_out.mean(axis=1).shape)
    #print(r_values.shape)
    # Plotear velocidad radial en el primer subplot
    axes[0].plot(r_values, vrad_out.mean(axis=1), '-', label=f'Snapshot {snapshot}', linewidth=1.5)
    axes[0].set_xlabel('Radial Distance', fontsize=14)
    axes[0].set_ylabel(r'Radial Velocity [ ]', fontsize=14)
    axes[0].set_title('Radial Velocity', fontsize=16)
    axes[0].legend(fontsize=10)
    axes[0].grid(True)

    # Plotear velocidad azimutal en el segundo subplot
    axes[1].plot(r_values, vazi_out.mean(axis=1), '-', label=f'Snapshot {snapshot}', linewidth=1.5)
    axes[1].set_xlabel('Radial Distance', fontsize=14)
    axes[1].set_ylabel(r'Azimuthal Velocity [ ]', fontsize=14)
    axes[1].set_title('Azimuthal Velocity', fontsize=16)
    axes[1].legend(fontsize=10)
    axes[1].grid(True)

plt.tight_layout()
plt.show()
