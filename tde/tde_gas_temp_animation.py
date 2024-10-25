import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import glob
import re

# Ruta donde está almacenado el archivo planet0.dat
path_sergei = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad_sergei/"
path_test = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad/"

path = path_sergei
outputs_para_plotear = 100
Ntot = outputs_para_plotear 


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

def cargar_temperature(output_number, path):
    dens    = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)*unit_density   #no units
    energy  = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX)*unit_energy    #no units
    press   = (gamma-1.0)*energy #no units
    temp    = mu*press/(dens*R_gas)# * unit_temperature
    return temp

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

#Ntot = 35 #len(valid_files)


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
vmin_gas, vmax_gas = None, None
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)
    log_dens_out = np.log10(dens_out)
    if vmin_gas is None or vmax_gas is None:
        vmin_gas = log_dens_out.min()
        vmax_gas = log_dens_out.max()
    else:
        vmin_gas = min(vmin_gas, log_dens_out.min())
        vmax_gas = max(vmax_gas, log_dens_out.max())


# Cargar la temperatura del frame 0 (T0)
T0 = cargar_temperature(0, path)
# Calcular los valores mínimos y máximos de la densidad globalmente
vmin_temp, vmax_temp = None, None
for output_number in range(Ntot):
    temp_out = cargar_temperature(output_number, path)
    log_temp_out = np.log10(temp_out)
    if vmin_temp is None or vmax_temp is None:
        vmin_temp = log_temp_out.min()
        vmax_temp = log_temp_out.max()
    else:
        vmin_temp = min(vmin_temp, log_temp_out.min())
        vmax_temp = max(vmax_temp, log_temp_out.max())

# Calcular los valores mínimos y máximos de la densidad globalmente
vmin_temp_diff, vmax_temp_diff = None, None
for output_number in range(Ntot):
    temp_out = cargar_temperature(output_number, path)
    Tdiff = (temp_out - T0)/T0
    log_temp_out_diff =Tdiff #np.log10(Tdiff)
    if vmin_temp_diff is None or vmax_temp_diff is None:
        vmin_temp_diff = log_temp_out_diff.min()
        vmax_temp_diff = log_temp_out_diff.max()
    else:
        vmin_temp_diff = min(vmin_temp_diff, log_temp_out_diff.min())
        vmax_temp_diff = max(vmax_temp, log_temp_out_diff.max())

print(vmin_temp_diff, vmax_temp_diff)

# Configuración de la figura con 3 subgráficos (densidad, temperatura y Tdiff)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 6))

# Gráfico de la densidad
c1 = ax1.pcolormesh(X, Y, np.zeros_like(X), cmap='jet', shading='nearest', vmin=vmin_gas, vmax=vmax_gas)
ax1.set_title('Density (log10)')
ax1.set_xlabel('X [AU]')
ax1.set_ylabel('Y [AU]')
cb1 = fig.colorbar(c1, ax=ax1, label=r'Density [$\mathrm{g/cm^2}$]')

# Gráfico de la temperatura
c2 = ax2.pcolormesh(X, Y, np.zeros_like(X), cmap='jet', shading='nearest', vmin=vmin_temp, vmax=vmax_temp)
ax2.set_title('Temperature (log10)')
ax2.set_xlabel('X [AU]')
ax2.set_ylabel('Y [AU]')
cb2 = fig.colorbar(c2, ax=ax2, label=r'Temperature [K]')

# Gráfico de la diferencia de temperatura Tdiff
c3 = ax3.pcolormesh(X, Y, np.zeros_like(X), cmap='bwr', shading='nearest', vmin=vmin_temp_diff, vmax=vmax_temp_diff)
ax3.set_title(r'Temperature Difference ($(T - T_0)/ T_0$)')
ax3.set_xlabel('X [AU]')
ax3.set_ylabel('Y [AU]')
cb3 = fig.colorbar(c3, ax=ax3, label=r'$(T - T_0)/ T_0$')

# Función para actualizar cada frame en la animación
def actualizar(frame):
    # Cargar la densidad y la temperatura para el frame actual
    dens_frame = cargar_densidad(frame, path)
    log_dens_frame = np.log10(dens_frame)
    temp_frame = cargar_temperature(frame, path)
    log_temp_frame = np.log10(temp_frame)
    
    # Calcular la diferencia de temperatura (Tdiff)
    Tdiff = (temp_frame - T0) / T0
    #print(Tdiff)

    # Actualizar los datos de densidad, temperatura y Tdiff
    c1.set_array(log_dens_frame.ravel())  # Actualiza los datos del pcolormesh de densidad
    c2.set_array(log_temp_frame.ravel())  # Actualiza los datos del pcolormesh de temperatura
    c3.set_array(Tdiff.ravel())           # Actualiza los datos del pcolormesh de Tdiff

    # Calcular el tiempo de evolución en años
    tiempo_evolucion = Ninter * DT * frame

    # Actualizar el título general con el tiempo de evolución
    fig.suptitle(f'Frame: {frame}, Evolution time: {tiempo_evolucion:.2e} yr', fontsize=16)

    # Guardar cada frame como PNG (opcional)
    file_path = os.path.join(output_dir, f"output_{frame:04d}.png")
    plt.savefig(file_path)

    return c1, c2, c3

# Crear la animación
ani = animation.FuncAnimation(fig, actualizar, frames=Ntot, interval=100, repeat=False)

# Guardar la animación en un archivo MP4 en el mismo directorio de los PNG
output_file = os.path.join(output_dir, 'gas_temp_animation_alfa1.mp4')
ani.save(output_file, writer='ffmpeg', fps=10)

# Mostrar la animación en pantalla
plt.show()