import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import glob
import re


#Datos de la simulacion
#path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d/"
#path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_iso/"
path = "/Users/matias/Simulations/mi_fargo3d/outputs/fargo/"
path = "/Users/matias/Simulations/mi_fargo3d/outputs/flyby2d_1MJ_YESfeel_lower/"
frame = 50



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

#unit_energy2 = unit_mass*(unit_length**2)/(unit_time*3.154e7)/unit_length**3    #erg/cm**3  #guido malo, falta tiempo2
#print(unit_energy, unit_energy2)

unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length)) #K


# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out*unit_surf_density

def cargar_temperature(output_number, path):
    dens    = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)*unit_density   #no units
    energy  = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX)*unit_energy    #no units
    press   = (gamma-1.0)*energy 
    
    #temp    = mu*press/(dens*R_gas) * unit_temperature
    temp = energy / dens / (R_gas / mu) * (gamma - 1)  # Temperatura en K
    
    return temp#*unit_temperature

# Función para leer las coordenadas del planeta desde planet0.dat
def leer_coordenadas_planeta(path, snapshot, file_planet):
    planet_file0 = path + file_planet
    planet0 = np.genfromtxt(planet_file0)
    xp = planet0[snapshot][1]  # Coordenada x del planeta
    yp = planet0[snapshot][2]  # Coordenada y del planeta
    return xp, yp

# Ruta donde está almacenado el archivo planet0.dat
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



time = frame * DT * Ninter
dens_frame = cargar_densidad(frame, path)
log_dens_frame = np.log10(dens_frame)

temp_frame = cargar_temperature(frame, path)
log_temp_frame = np.log10(temp_frame)

# Crear la figura y las subparcelas
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 8))

# Gráfico de la densidad
c1 = ax1.pcolormesh(X, Y, log_dens_frame, cmap='jet', shading='nearest', vmin=vmin_gas, vmax=vmax_gas)
ax1.set_title('Density (log10)')
ax1.set_xlabel('X [AU]')
ax1.set_ylabel('Y [AU]')
fig.colorbar(c1, ax=ax1, label=r'Density [$\mathrm{g/cm^2}$]')
# Gráfico de la temperatura
c2 = ax2.pcolormesh(X, Y, log_temp_frame, cmap='jet', shading='nearest', vmin=vmin_temp, vmax=vmax_temp)
ax2.set_title('Temperature (log10)')
ax2.set_xlabel('X [K]')
ax2.set_ylabel('Y [K]')
fig.colorbar(c2, ax=ax2, label='Temperature [K]')

# Agregar el título general con el tiempo de evolución
time_label = f'Evolution time: {time:.2f} años'  # Reemplaza 'time' con la variable correspondiente al tiempo
fig.suptitle(time_label, fontsize=16)

# Mostrar la figura
plt.tight_layout()
plt.show()

