import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import glob
import re

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out

# Función para leer las coordenadas del planeta desde planet0.dat
def leer_coordenadas_planeta(path, snapshot, file_planet):
    planet_file0 = path + file_planet
    planet0 = np.genfromtxt(planet_file0)
    xp = planet0[snapshot][1]  # Coordenada x del planeta
    yp = planet0[snapshot][2]  # Coordenada y del planeta
    return xp, yp

# Ruta donde está almacenado el archivo planet0.dat
path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d/"
planet_file = path + "planet0.dat"
planet_data = np.genfromtxt(planet_file)

# Extraer las coordenadas x e y del planeta
xp = planet_data[:, 1]  # Coordenada x del planeta
yp = planet_data[:, 2]  # Coordenada y del planeta

# Número de snapshots disponibles
variables_par = np.genfromtxt(path+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo var    iable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                         #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                              #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                                  #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                     #

def P(parametro):
        return valores_par[parametros_par.index(parametro)] 

# Identificar el número de outputs disponibles
# Patrón para archivos gasdens*.dat con un número entero
pattern = re.compile(r'gasdens(\d+)\.dat')
# Lista de archivos que coinciden con el patrón gasdens*.dat
files = glob.glob(path + "gasdens*.dat")
# Filtrar archivos que se ajustan al patrón correcto
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
# Contar el número de archivos válidos
Ntot = len(valid_files)

Ninter = int(P("NINTERM"))
DT = float(P("DT"))

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = int(P("NX")) #len(domain_x) - 1
NY = int(P("NY")) # len(domain_y) - 1


#cargar datos:
# Calcular los valores mínimos y máximos de la densidad globalmente
vmin, vmax = None, None
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)
    log_dens_out = np.log10(dens_out)
    if vmin is None or vmax is None:
        vmin = log_dens_out.min()
        vmax = log_dens_out.max()
    else:
        vmin = min(vmin, log_dens_out.min())
        vmax = max(vmax, log_dens_out.max())

# Configuración de la grilla
def Grilla_XY():
    R = 0.5 * (domain_y[1:] + domain_y[:-1])
    Phi = 0.5 * (domain_x[1:] + domain_x[:-1])
    # Aseguramos que la grilla incluye el valor final para cerrar el círculo
    #Phi = np.append(Phi, Phi[0] + 2 * np.pi)  # Añadir el primer valor + 2pi al final #ojo que con esto no me permite graficar pcolomersh usando estos X e Y
    P, R = np.meshgrid(Phi, R)
    X, Y = R * np.cos(P), R * np.sin(P)
    return X, Y

X, Y = Grilla_XY()

# Función para calcular el radio de Hill
def calcular_radio_hill(r_planeta, masa_planeta, masa_estrella):
    return r_planeta * (masa_planeta / (3 * masa_estrella))**(1/3)


# Parámetros del planeta y la estrella
masa_planeta = 1e-3  # Masa del planeta en unidades de masas solares (por ejemplo, Tierra en comparación con el Sol)
masa_estrella = 1.0  # Masa de la estrella en unidades solares (por ejemplo, 1 para el Sol)
r_planeta = np.sqrt(xp**2 + yp**2)  # Calcular el radio del planeta en coordenadas polares

# Calcular el radio de Hill
radio_hill = calcular_radio_hill(r_planeta, masa_planeta, masa_estrella)


frame = 2
dens_frame = cargar_densidad(frame, path)
log_dens_out = np.log10(dens_frame)


# Función para plotear la grilla
def plot_grilla(X, Y):
    plt.figure(figsize=(8, 8))
    plt.plot(X, Y, color='gray', lw=0.5)  # Lineas horizontales (radiales)
    plt.plot(X.T, Y.T, color='gray', lw=0.5)  # Lineas verticales (azimutales)
    plt.scatter(X, Y, color='red', s=1)  # Puntos de la grilla
    #plt.scatter(xp, yp, color='blue', s=100, label="Planeta")  # Posición del planeta
    circle = plt.Circle((xp[0], yp[0]), radio_hill[0], color='black', fill=False, linestyle='--', label="Radio de Hill")
    plt.gca().add_patch(circle)
    plt.pcolormesh(X, Y, log_dens_out, cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
    plt.xlabel('X (AU)')
    plt.ylabel('Y (AU)')
    plt.title('Grilla en coordenadas cartesianas (X, Y)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# Llamamos a la función para plotear la grilla
plot_grilla(X, Y)


frame = Ntot - 1
dens_out = cargar_densidad(frame, path)
log_dens_out = np.log10(dens_out)

fig = plt.figure(figsize=(8, 8))
plt.pcolormesh(X, Y, log_dens_out, cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
plt.show()