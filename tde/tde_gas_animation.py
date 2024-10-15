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
x_coords = planet_data[:, 1]  # Coordenada x del planeta
y_coords = planet_data[:, 2]  # Coordenada y del planeta

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
Rmax = float(P("YMAX"))

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = len(domain_x) - 1
NY = len(domain_y) - 1

# Configuración de la grilla
def Grilla_XY():
    R = 0.5 * (domain_y[1:] + domain_y[:-1])
    Phi = 0.5 * (domain_x[1:] + domain_x[:-1])
    P, R = np.meshgrid(Phi, R)
    X, Y = R * np.cos(P), R * np.sin(P)
    return X, Y

X, Y = Grilla_XY()

# Crear la carpeta para guardar los PNGs si no existe
output_dir = os.path.join(path, "gas_png")
os.makedirs(output_dir, exist_ok=True)

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

# Función para actualizar cada frame en la animación y guardar PNGs
def actualizar(frame):
    plt.gca().cla()  # Limpiar el gráfico para el nuevo frame
    
    # Cargar la densidad de gas para el frame actual
    dens_out = cargar_densidad(frame, path)
    log_dens_out = np.log10(dens_out)

    # Graficar la densidad de gas
    plt.pcolormesh(X, Y, log_dens_out, cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
    
    # Graficar la órbita completa del planeta
    #plt.plot(x_coords, y_coords, color='black', linestyle='-', marker='o', markersize=5, label='Órbita del planeta')
    
    # Graficar la posición actual del planeta en el frame
    #plt.scatter(x_coords[frame], y_coords[frame], color='teal', s=100, edgecolor='black', label=f'Posición en frame {frame+1}')
    
    # Graficar la posición inicial en azul
    #plt.scatter(x_coords[0], y_coords[0], color='blue', s=100, label='Posición inicial')
    
    # Graficar la posición final en rojo
    #plt.scatter(x_coords[-1], y_coords[-1], color='red', s=100, label='Posición final')
    
    # Graficar el Sol en el centro
    #plt.scatter(0, 0, color='yellow', s=200, edgecolor='orange', label='Star', zorder=5)

    # Calcular el tiempo de evolución en años
    tiempo_evolucion = Ninter * DT * frame

    # Configurar etiquetas, límites y leyenda
    plt.xlim([-Rmax, Rmax])
    plt.ylim([-Rmax, Rmax])
    plt.xlabel('X [AU]', fontsize=14)
    plt.ylabel('Y [AU]', fontsize=14)
    plt.title(f'Frame {frame} \nTime {tiempo_evolucion:.2e} yr', fontsize=16)

    plt.legend()
    plt.gca().set_aspect('equal', adjustable='box')  # Mantener proporciones
    plt.grid(False)

    # Guardar cada frame como PNG
    file_path = os.path.join(output_dir, f"output_{frame:04d}.png")
    plt.savefig(file_path)

# Crear la animación
fig = plt.figure(figsize=(8, 8))
ani = animation.FuncAnimation(fig, actualizar, frames=Ntot, interval=100, repeat=False)

# Guardar la animación en un archivo MP4 en el mismo directorio de los PNG
output_file = os.path.join(output_dir, 'planet_orbit_with_gas_animation.mp4')
ani.save(output_file, writer='ffmpeg', fps=10)

# Mostrar la animación en pantalla
plt.show()
