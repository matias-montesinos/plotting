import numpy as np
import matplotlib.pyplot as plt
import os
import sys



# Ruta donde está almacenado el archivo planet0.dat
path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad/"

variables_par = np.genfromtxt(path+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo var    iable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                         #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                              #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                                  #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                     #

def P(parametro):
        return valores_par[parametros_par.index(parametro)] 


#UNITS
gamma=1.6666
mu=2.35
kb=1.380650424e-16 #ergs
R_gas=8.314472e7 #erg/K/mol 
G = 6.67259e-8 #cm**3 g**-1 s**-2
mp= 1.672623099e-24 # g
unit_mass   = 1.9891e33 # g
unit_length = 1.49598e13 #cm

unit_density = (unit_mass)/(unit_length**3) # g/cm**3
unit_surf_density = (unit_mass/unit_length**2) # g/cm**2
unit_time = np.sqrt( pow(unit_length,3.) / G / unit_mass)/ 3.154e7 #yr 

unit_energy = unit_mass/(unit_time*3.154e7)**2/(unit_length) #erg/cm3 = gr/s2/cm
unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length)) #K

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out*unit_surf_density

# Función para cargar la temperatura de gas
def cargar_temperature(output_number, path):
    dens = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)*unit_density
    energy = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX)*unit_energy
    press = (gamma-1.0)*energy
    temp = mu*press/(dens*R_gas)
    return temp

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]  # Mantener los límites fantasma

NX = int(P("NX")) #len(domain_x) - 1
NY = int(P("NY")) #len(domain_y) - 1
Ninter = int(P("NINTERM"))
DT = float(P("DT"))

# Solicitar al usuario el snapshot a visualizar
snapshot = int(input("Ingrese el número del snapshot: "))

tiempo_evolucion = Ninter * DT * snapshot




# Cargar los datos de densidad y temperatura del snapshot seleccionado
densidad = cargar_densidad(snapshot, path)
temperatura = cargar_temperature(snapshot, path)

# Tomar el promedio azimutal para obtener perfiles radiales
densidad_promedio = np.mean(densidad, axis=1)
temperatura_promedio = np.mean(temperatura, axis=1)

# Asegurarse de que la grilla y los datos tengan la misma longitud
domain_y = domain_y[:NY]

# Graficar los perfiles radiales de densidad y temperatura
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot de densidad con un color más moderno y línea más gruesa con transparencia
ax1.plot(domain_y, densidad_promedio, color='#2b6ca3', linestyle='-', label="Density", lw=3, alpha=0.8)
ax1.set_xlabel('Radius [AU]', fontsize=16)
ax1.set_ylabel('Density [$g/cm^2$]', color='#2b6ca3', fontsize=16)

# Ajustar el tamaño de las etiquetas de los ticks en el eje X y Y
ax1.tick_params(axis='both', labelsize=14)
ax1.tick_params(axis='y', labelcolor='#2b6ca3')

# Crear un segundo eje y para la temperatura
ax2 = ax1.twinx()

# Cambiar el color de la temperatura a un color más agradable (verde suave)
ax2.plot(domain_y, temperatura_promedio, color='#d35400', linestyle='-', label="Temperature", lw=3, alpha=0.8)
ax2.set_ylabel('Temperature [K]', color='#d35400', fontsize=16)

# Ajustar el tamaño de las etiquetas de los ticks en el eje Y
ax2.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='y', labelcolor='#d35400')

fig.suptitle(f'Frame: {snapshot}, Evolution time: {tiempo_evolucion:.2e} yr', fontsize=16)

# Mostrar la gráfica
plt.tight_layout()
plt.show()
