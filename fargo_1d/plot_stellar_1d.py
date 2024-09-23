import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
from PIL import Image
import glob
import re


path = "/home/matias/Simulations/fargo3d/outputs/accretion/"

variables_par = np.genfromtxt(path+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo var    iable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                         #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                              #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                                  #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                     #

def P(parametro):
        return valores_par[parametros_par.index(parametro)]                                                                                                                                                                                 #Retorna siempre str, recordad tranformar si necesita un int o float


# Identificar el número de outputs disponibles
# Patrón para archivos gasdens*.dat con un número entero
pattern = re.compile(r'gasdens(\d+)\.dat')
# Lista de archivos que coinciden con el patrón gasdens*.dat
files = glob.glob(path + "gasdens*.dat")
# Filtrar archivos que se ajustan al patrón correcto
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
# Contar el número de archivos válidos
Ntot = len(valid_files)
print(Ntot)


# Datos de la grilla
Rmin = float(P("YMIN"))
Rmax = float(P("YMAX"))

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = len(domain_x) - 1
NY = len(domain_y) - 1
dr = (Rmax - Rmin) / NY

# Coordenada r
r_values = np.arange(Rmin, Rmax, dr)

Rstar = Rmin
u_func = -Rstar/r_values
#u_func = r_values/Rstar

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

# Solicitar al usuario el snapshot a usar
snapshot = 1 # int(input("Ingrese el número de snapshot: "))
dens_out = cargar_densidad(snapshot, path)
vrad_out = cargar_vrad(snapshot, path)
vazi_out = cargar_vazi(snapshot, path)
print(dens_out)

# Configuración de los gráficos
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

# Primer subplot para vrad_out
axes[0].plot(u_func, vrad_out, '-', label='radial velocity', linewidth=2)
axes[0].set_xlabel('Radial Distance', fontsize=14)
axes[0].set_ylabel(r'Radial Velocity [ ]', fontsize=14)
axes[0].set_title('Radial Velocity', fontsize=16)
axes[0].legend(fontsize=12)
axes[0].grid(True)

# Segundo subplot para vazi_out
axes[1].plot(r_values, vazi_out, '-', label='azimuthal velocity', linewidth=2)
axes[1].set_xlabel('Radial Distance', fontsize=14)
axes[1].set_ylabel(r'Azimuthal Velocity [ ]', fontsize=14)
axes[1].set_title('Azimuthal Velocity', fontsize=16)
axes[1].legend(fontsize=12)
axes[1].grid(True)

plt.tight_layout()
plt.show()