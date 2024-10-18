import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import glob
import re


# Ruta donde está almacenado el archivo planet0.dat
path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad/"

#UNITS
gamma=1.6666
mu=2.35 #
kb=1.380650424e-16 #ergs
R_gas=8.314472e7 #erg/K/mol 
G = 6.67259e-8 #cm**3 g**-1 s**-2
mp= 1.672623099e-24 # g
unit_mass   = 1.9891e33 # g
unit_length = 1.49598e13 #cm

unit_density=(unit_mass)/(unit_length**3) # g/cm**3
unit_surf_density = (unit_mass/unit_length**2) # g/cm**2
unit_time     = np.sqrt( pow(unit_length,3.) / G / unit_mass)/ 3.154e7 #yr 

unit_energy   = unit_mass/(unit_time*3.154e7)**2/(unit_length)                  #erg/cm3 = gr/s2/cm  #juan garrido
unit_surf_energy = unit_energy / unit_length

unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length)) #K

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out*unit_surf_density

# Función para leer las coordenadas del planeta desde planet0.dat
def leer_coordenadas_planeta(path, snapshot, file_planet):
    planet_file0 = path + file_planet
    planet0 = np.genfromtxt(planet_file0)
    xp = planet0[snapshot][1]  # Coordenada x del planeta
    yp = planet0[snapshot][2]  # Coordenada y del planeta
    return xp, yp



planet_file = path + "planet0.dat"
planet_data = np.genfromtxt(planet_file)

# Extraer las coordenadas x e y del planeta
x_coords = planet_data[:, 1]  # Coordenada x del planeta
y_coords = planet_data[:, 2]  # Coordenada y del planeta


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
# Crear la animación
fig, ax = plt.subplots(figsize=(8, 8))

# Graficar la densidad de gas en el primer frame (solo para configurar la colorbar)
dens_out = cargar_densidad(0, path)
log_dens_out = np.log10(dens_out)
plot = ax.pcolormesh(X, Y, log_dens_out, cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)

# Agregar la colorbar una sola vez, fuera del loop de actualización
cbar = fig.colorbar(plot, ax=ax)
cbar.set_label(r'Log Density [$\mathrm{g/cm^2}$]', fontsize=14)  # Etiqueta de la colorbar

# Función para actualizar cada frame en la animación
def actualizar(frame):
    ax.cla()  # Limpiar el gráfico para el nuevo frame
    
    # Cargar la densidad de gas para el frame actual
    dens_out = cargar_densidad(frame, path)
    log_dens_out = np.log10(dens_out)

    # Graficar la densidad de gas
    plot = ax.pcolormesh(X, Y, log_dens_out, cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
    
    # Calcular el tiempo de evolución en años
    tiempo_evolucion = Ninter * DT * frame

    # Configurar etiquetas, límites y leyenda
    ax.set_xlim([-Rmax, Rmax])
    ax.set_ylim([-Rmax, Rmax])
    ax.set_xlabel('X [AU]', fontsize=14)
    ax.set_ylabel('Y [AU]', fontsize=14)
    ax.set_title(f'Frame {frame} \nTime {tiempo_evolucion:.2e} yr', fontsize=16)

    ax.set_aspect('equal', adjustable='box')  # Mantener proporciones
    ax.grid(False)

    # Guardar cada frame como PNG
    file_path = os.path.join(output_dir, f"output_{frame:04d}.png")
    plt.savefig(file_path)

# Crear la animación
ani = animation.FuncAnimation(fig, actualizar, frames=Ntot, interval=100, repeat=False)

# Guardar la animación en un archivo MP4 en el mismo directorio de los PNG
output_file = os.path.join(output_dir, 'gas_isothermal.mp4')
ani.save(output_file, writer='ffmpeg', fps=10)

# Mostrar la animación en pantalla
plt.show()
