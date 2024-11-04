import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
import glob
import re


path_tde_sergeiElite = "/home/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad_sergei/" # frame 100
path_tde_segei_imac = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad_sergei/"
path = path_tde_segei_imac

# Identificar el número de outputs disponibles
# Patrón para archivos gasdens*.dat con un número entero
pattern = re.compile(r'gasdens(\d+)\.dat')
# Lista de archivos que coinciden con el patrón gasdens*.dat
files = glob.glob(path + "gasdens*.dat")
# Filtrar archivos que se ajustan al patrón correcto
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
# Contar el número de archivos válidos
Ntot = len(valid_files)
print("En total hay: ", Ntot, "outputs")

#UNITS
gamma=1.6666
mu=2.35 #
kb=1.380650424e-16 #ergs/K
R_gas=8.314472e7 #erg/K/mol 
G = 6.67259e-8 #cm**3 g**-1 s**-2
mp= 1.672623099e-24 # g

unit_mass   = 1.9891e33 # g
unit_length = 1.49598e13 #cm
unit_density=(unit_mass)/(unit_length**3) # g/cm**3
unit_surf_density = (unit_mass/unit_length**2) # g/cm**2
unit_time     = np.sqrt( pow(unit_length,3.) / G / unit_mass) #s 
unit_energy   = unit_mass/(unit_time)**2/(unit_length)  #erg/cm3 = gr/s2/cm  #juan garrido
unit_surf_energy = unit_energy / unit_length
unit_velocity = unit_length / unit_time  # cm/s

unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length)) #K


variables_par = np.genfromtxt(path+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo var    iable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                         #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                              #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                                  #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                     #

def P(parametro):
        return valores_par[parametros_par.index(parametro)] 



def compute_accretion_rate(output_number, path, r0):
    """
    Calcula la tasa de acreción en el radio r0 para un output_number dado.

    Parámetros:
    - output_number: número del output (snapshot) de la simulación.
    - path: ruta al directorio donde están los datos de la simulación.
    - r0: radio seleccionado donde calcular la tasa de acreción (en unidades de AU).

    Retorna:
    - accretion_rate: tasa de acreción en unidades de masas solares por año (Msun/yr).
    """

    # Cargar los datos de densidad superficial (gasdens) y velocidad radial (gasvy)
    dens_data = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX) * unit_surf_density  # [g/cm^2]
    vr_data = np.fromfile(path + f"gasvy{output_number}.dat").reshape(NY, NX) * unit_velocity  # [cm/s]

    # Obtener el radio y ángulo azimutal en cada celda
    Phi = 0.5 * (domain_x[1:] + domain_x[:-1])  # [rad]
    R = 0.5 * (domain_y[1:] + domain_y[:-1])    # [AU]

    # Si la malla es uniforme en phi
    Delta_phi = 2 * np.pi / NX  # [rad]

    # Encontrar el índice radial más cercano a r0
    radial_index = np.argmin(np.abs(R - r0))
    r_selected = R[radial_index] * unit_length  # Convertir a [cm]

    # Extraer densidad y velocidad radial en el radio seleccionado
    Sigma_r = dens_data[radial_index, :]  # [g/cm^2]
    vr_r = vr_data[radial_index, :]       # [cm/s]

    # Calcular el flujo de masa en cada celda azimutal
    mass_flux = -Sigma_r * vr_r * r_selected * Delta_phi  # [g/s]

    # Sumar sobre todas las celdas azimutales para obtener la tasa de acreción total
    accretion_rate = np.sum(mass_flux)  # [g/s]

    # Convertir la tasa de acreción a unidades de masas solares por año
    accretion_rate_msun_per_year = accretion_rate * (3.154e7) / 1.989e33  # [Msun/yr]

    return accretion_rate_msun_per_year




# construyendo la grilla
Ninter = int(P("NINTERM"))
DT = float(P("DT"))

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = int(P("NX"))# len(domain_x) - 1
NY = int(P("NY")) #len(domain_y) - 1


output_number = 1
r0 = 1
mdot = compute_accretion_rate(output_number, path, r0)
output_numbers = range(Ntot)
print(f"Tasa de acreción en r = {r0} AU: {mdot:.3e} M☉/año")


# Lista para almacenar la tasa de acreción y los tiempos
mdot_list = []
time_list = []

# Calcula la tasa de acreción para cada output
for output_number in output_numbers:
    # Tiempo en unidades de código
    time_code_units = output_number * Ninter * DT
    # Convertir el tiempo a segundos y luego a días
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    # Calcular la tasa de acreción
    mdot = compute_accretion_rate(output_number, path, r0)
    mdot_list.append(mdot)
    time_list.append(time_days)

# Convertir listas a arrays para facilitar el manejo
mdot_array = np.array(mdot_list)
time_array = np.array(time_list)

# Grafica la tasa de acreción en función del tiempo (en días)
plt.figure(figsize=(10, 6))
plt.plot(time_array, mdot_array, marker='o')
plt.xlabel('Time [days]')
plt.ylabel('Accretion rate [M$_\odot$/año]')
plt.title(f'Accretion rate at r = {r0} AU')
plt.grid(True)
plt.show()


# Define the range of radii
# Create an array of radii
radii = 0.5 * (domain_y[1:] + domain_y[:-1]) 
num_radii = NY

# Assuming 'output_numbers' and 'time_array' are already defined
# output_numbers = range(Ntot)
# time_array = [...]  # Calculated or existing time array in days

# Initialize a matrix to store \dot{M}(r, t)
mdot_matrix = np.zeros((len(time_array), num_radii))  # Dimensions: time x radius

# Loop over times (outputs)
for i, output_number in enumerate(output_numbers):
    # Optional: print progress
    print(f"Processing output {output_number + 1} of {Ntot}")
    
    # Loop over radii
    for j, r0 in enumerate(radii):
        # Calculate the accretion rate at radius r0 and time t
        mdot = compute_accretion_rate(output_number, path, r0)
        # Store in the matrix
        mdot_matrix[i, j] = mdot


# Generate the plot of \dot{M}(r, t) as a heatmap
# Create a meshgrid of radii and times for plotting
R, T = np.meshgrid(radii, time_array)  # R in [AU], T in [days]

# Create the figure and axis
plt.figure(figsize=(12, 6))

# Use pcolormesh for the heatmap
plt.pcolormesh(R, T, mdot_matrix, shading='gouraud', cmap='viridis')

# Add a color bar
cbar = plt.colorbar()
cbar.set_label('Accretion rate [M$_\odot$/yr]')  # Color bar label in English

# Labels and title in English
plt.xlabel('Radius [AU]')
plt.ylabel('Time [days]')
plt.title('Accretion rate $\dot{M}(r, t)$')

# Optionally adjust the limits of the plot
plt.xlim(radii.min(), radii.max())
plt.ylim(time_array.min(), time_array.max())

# Show the plot
plt.tight_layout()
plt.show()


############
# plot de Mdot(r) para tistintos tiempos
###########
# Seleccionar los tiempos (outputs) en los que deseas graficar \dot{M}(r)
# Por ejemplo, puedes seleccionar 3 tiempos: al inicio, a la mitad y al final de la simulación
selected_outputs = [0, 10, 50, 100]  # Puedes ajustar estos valores

# Crear una lista con los tiempos correspondientes en días
selected_times = []
for output_number in selected_outputs:
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    selected_times.append(time_days)

# Crear un array de radios
radii =  0.5 * (domain_y[1:] + domain_y[:-1])  # [AU] 
# Inicializar una lista para almacenar las curvas de \dot{M}(r) en los tiempos seleccionados
mdot_curves = []

# Bucle sobre los outputs seleccionados
for idx, output_number in enumerate(selected_outputs):
    # Opcional: imprimir el progreso
    print(f"Processing output {output_number + 1} of {Ntot} for time {selected_times[idx]:.1f} days")
    
    # Inicializar un array para \dot{M}(r) en este tiempo
    mdot_r = []
    # Bucle sobre los radios
    for r0 in radii:
        # Calcular la tasa de acreción en el radio r0 y tiempo t
        mdot = compute_accretion_rate(output_number, path, r0)
        mdot_r.append(mdot)
    # Agregar la curva a la lista
    mdot_curves.append(mdot_r)

# Graficar las curvas de \dot{M}(r) para los tiempos seleccionados
plt.figure(figsize=(10, 6))

for i, mdot_r in enumerate(mdot_curves):
    plt.plot(radii, mdot_r, label=f'Time = {selected_times[i]:.1f} days')

plt.xlabel('Radius [AU]')
plt.ylabel('Accretion rate [M$_\odot$/yr]')
plt.title('Accretion rate $\dot{M}(r)$ at different times')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


## velocidades radiales
# Crear una lista con los tiempos correspondientes en días
selected_times = []
for output_number in selected_outputs:
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    selected_times.append(time_days)

# Cargar la grilla si no lo has hecho
#domain_x = np.genfromtxt(path + "domain_x.dat")
#domain_y = np.genfromtxt(path + "domain_y.dat")
#NX = len(domain_x) - 1
#NY = len(domain_y) - 1

# Calcular los centros de las celdas radiales
radii = 0.5 * (domain_y[1:] + domain_y[:-1])  # [AU]

# Inicializar una lista para almacenar las velocidades radiales promediadas
vr_profiles = []

# Bucle sobre los outputs seleccionados
for idx, output_number in enumerate(selected_outputs):
    print(f"Procesando output {output_number} para tiempo {selected_times[idx]:.1f} días")
    
    # Cargar los datos de velocidad radial
    vr_data = np.fromfile(path + f"gasvy{output_number}.dat").reshape(NY, NX) * unit_velocity  # [cm/s]
    
    # Promediar la velocidad radial en la dirección azimutal
    vr_mean = np.mean(vr_data, axis=1)  # Promedio sobre el eje azimutal (axis=1)
    
    # Almacenar el perfil promediado
    vr_profiles.append(vr_mean)

# Convertir la velocidad a km/s si lo deseas
vr_profiles_km_s = [vr / 1e5 for vr in vr_profiles]  # Convertir de cm/s a km/s

# Graficar los perfiles de velocidad radial para los tiempos seleccionados
plt.figure(figsize=(10, 6))

for i, vr_mean in enumerate(vr_profiles_km_s):
    plt.plot(radii, vr_mean, label=f'Time = {selected_times[i]:.1f} days')

plt.xlabel('Radius [AU]')
plt.ylabel('Radial Velocity [km/s]')
plt.title('Azimuthally Averaged Radial Velocity at Different Times')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



##
### Calculo de la luminosidad de acrecion

# Parámetros del objeto central y constantes
M_star = 1.0  # Masa del objeto central en masas solares
M_star_cgs = M_star * 1.989e33  # Convertir a gramos
G_cgs = 6.67430e-8  # Constante de gravitación universal en cm^3 g^-1 s^-2
L_sun_cgs = 3.828e33  # Luminosidad solar en erg/s
M_sun = 4.83  # Magnitud absoluta del Sol en la banda V

# Convertir radios a cm
radii_cm = radii * unit_length  # unit_length en cm

# Cálculo de la luminosidad total y su conversión a magnitudes
time_array_L = []
L_total_Lsun = []
M_total = []

for idx, output_number in enumerate(output_numbers):
    # Calcular el tiempo correspondiente
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    time_array_L.append(time_days)
    
    # Calcular \dot{M}(r, t)
    mdot_r = []
    for r0 in radii:
        mdot = compute_accretion_rate(output_number, path, r0)
        mdot_r.append(mdot)
    mdot_r = np.array(mdot_r)
    
    # Convertir \dot{M} a g/s
    mdot_r_cgs = mdot_r * (1.989e33) / (3.154e7)  # g/s
    
    # Establecer a cero los valores negativos de \dot{M}
    mdot_r_cgs = np.where(mdot_r_cgs > 0, mdot_r_cgs, 0.0)
    
    # Verificar si \dot{M} es cero en todos los radios
    if np.all(mdot_r_cgs == 0):
        print(f"At time {time_days:.1f} days, \dot{{M}} is zero or negative at all radii.")
        L_total_Lsun_value = 0.0
        M = np.nan  # O asignar un valor muy grande
    else:
        # Calcular dL/dr
        dLdr = (3 * G_cgs * M_star_cgs * mdot_r_cgs) / (2 * (radii_cm ** 2))  # erg/s/cm
        
        # Integrar para obtener L_total
        L_total = np.trapz(dLdr, radii_cm)  # erg/s
        L_total_Lsun_value = L_total / L_sun_cgs  # Convertir a L_sun
        
        # Calcular la magnitud absoluta
        M = M_sun - 2.5 * np.log10(L_total_Lsun_value) + 17.0
    
    # Almacenar los valores
    L_total_Lsun.append(L_total_Lsun_value)
    M_total.append(M)

# Convertir listas a arrays numpy
time_array_L = np.array(time_array_L)
L_total_Lsun = np.array(L_total_Lsun)
M_total = np.array(M_total)

# Graficar la magnitud absoluta vs. tiempo
plt.figure(figsize=(10, 6))
plt.plot(time_array_L, L_total_Lsun, marker='o')
plt.xlabel('Time [days]')
#plt.ylabel('Absolute Magnitude $M$')
plt.ylabel('Accretion luminosity (bolometric) $L_\odot$')
plt.title('Total Accretion Luminosity Magnitude $M(t)$')
#plt.gca().invert_yaxis()  # Invertir eje Y
plt.grid(True)
plt.tight_layout()
plt.show()



### calculo de la luminosidad con extincion Av = 17 mag
# Definir la extinción en la banda V
# Definir la extinción en la banda V
A_V = 17.0  # magnitudes

# Parámetros del objeto central y constantes
M_sun = 4.83  # Magnitud absoluta del Sol en la banda V

# Inicializar listas para almacenar los resultados
time_array_L = []
L_total_Lsun = []
M_total = []

for idx, output_number in enumerate(output_numbers):
    # Calcular el tiempo correspondiente
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    time_array_L.append(time_days)
    
    # Calcular \dot{M}(r, t)
    mdot_r = []
    for r0 in radii:
        mdot = compute_accretion_rate(output_number, path, r0)
        mdot_r.append(mdot)
    mdot_r = np.array(mdot_r)
    
    # Convertir \dot{M} a g/s
    mdot_r_cgs = mdot_r * (1.989e33) / (3.154e7)  # g/s
    
    # Establecer a cero los valores negativos de \dot{M}
    mdot_r_cgs = np.where(mdot_r_cgs > 0, mdot_r_cgs, 0.0)
    
    # Verificar si \dot{M} es cero en todos los radios
    if np.all(mdot_r_cgs == 0):
        print(f"At time {time_days:.1f} days, \dot{{M}} is zero or negative at all radii.")
        L_total_Lsun_value = 0.0
        M_obs = np.nan  # Definir M_obs como NaN
    else:
        # Calcular dL/dr
        dLdr = (3 * G_cgs * M_star_cgs * mdot_r_cgs) / (2 * (radii_cm ** 2))  # erg/s/cm
        
        # Integrar para obtener L_total
        L_total = np.trapz(dLdr, radii_cm)  # erg/s
        L_total_Lsun_value = L_total / L_sun_cgs  # Convertir a L_sun
        
        # Calcular la magnitud absoluta intrínseca
        M_intrinsec = M_sun - 2.5 * np.log10(L_total_Lsun_value)
        
        # Calcular la magnitud observada con extinción
        M_obs = M_intrinsec + A_V
    
    # Almacenar los valores
    L_total_Lsun.append(L_total_Lsun_value)
    M_total.append(M_obs)


# Convertir listas a arrays numpy
time_array_L = np.array(time_array_L)
L_total_Lsun = np.array(L_total_Lsun)
M_total = np.array(M_total)

# Graficar la magnitud absoluta total vs. tiempo con extinción
plt.figure(figsize=(10, 6))
plt.plot(time_array_L, M_total, marker='o', label=f'Extinción $A_V={A_V}$ mag')
plt.xlabel('Time [days]')
plt.ylabel('Absolute Magnitude $M$')
plt.title('Total Accretion Luminosity Magnitude $M(t)$ with Extinction $A_V=17$ mag')
plt.gca().invert_yaxis()  # Invertir eje Y para que magnitudes más bajas (más brillantes) estén arriba
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



### ahora hacemos una correccion de la magnitud a la banda K
##

# Definir la extinción en la banda V
A_V = 17.0  # magnitudes

# Calcular la extinción en la banda K usando una ley de extinción estándar
A_K = 0.112 * A_V  # magnitudes

# Magnitud absoluta del Sol en la banda K
M_sun_K = 3.28  # magnitud absoluta del Sol en la banda K

# Inicializar listas para almacenar los resultados
time_array_L = []
L_total_Lsun = []
M_total_K = []

for idx, output_number in enumerate(output_numbers):
    # Calcular el tiempo correspondiente
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    time_array_L.append(time_days)
    
    # Calcular \dot{M}(r, t)
    mdot_r = []
    for r0 in radii:
        mdot = compute_accretion_rate(output_number, path, r0)
        mdot_r.append(mdot)
    mdot_r = np.array(mdot_r)
    
    # Convertir \dot{M} a g/s
    mdot_r_cgs = mdot_r * (1.989e33) / (3.154e7)  # g/s
    
    # Establecer a cero los valores negativos de \dot{M}
    mdot_r_cgs = np.where(mdot_r_cgs > 0, mdot_r_cgs, 0.0)
    
    # Verificar si \dot{M} es cero en todos los radios
    if np.all(mdot_r_cgs == 0):
        print(f"At time {time_days:.1f} days, \dot{{M}} is zero or negative at all radii.")
        L_total_Lsun_value = 0.0
        M_obs_K = np.nan  # Magnitud indefinida o extremadamente débil
    else:
        # Calcular dL/dr
        dLdr = (3 * G_cgs * M_star_cgs * mdot_r_cgs) / (2 * (radii_cm ** 2))  # erg/s/cm
        
        # Integrar para obtener L_total
        L_total = np.trapz(dLdr, radii_cm)  # erg/s
        L_total_Lsun_value = L_total / L_sun_cgs  # Convertir a L_sun
        
        # Calcular la magnitud absoluta intrínseca en la banda K
        M_intrinsec_K = M_sun_K - 2.5 * np.log10(L_total_Lsun_value)
        
        # Calcular la magnitud observada con extinción en la banda K
        M_obs_K = M_intrinsec_K + A_K
    
    # Almacenar los valores
    L_total_Lsun.append(L_total_Lsun_value)
    M_total_K.append(M_obs_K)

# Convertir listas a arrays numpy
time_array_L = np.array(time_array_L)
L_total_Lsun = np.array(L_total_Lsun)
M_total_K = np.array(M_total_K)

# Graficar la magnitud absoluta total vs. tiempo con extinción en la banda K
plt.figure(figsize=(10, 6))
plt.plot(time_array_L, M_total_K, marker='o', label=f'Extinción $A_K={A_K:.2f}$ mag')
plt.xlabel('Time [days]')
plt.ylabel('Absolute Magnitude $M_K$')
plt.title('Total Accretion Luminosity Magnitude $M_K(t)$ con Extinción $A_K=1.90$ mag')
plt.gca().invert_yaxis()  # Invertir eje Y para que magnitudes más bajas (más brillantes) estén arriba
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


## calculo de la luminosidad considerando extincion
# Definir la extinción en la banda V
A_V = 12.0  # magnitudes

# Definir la función para aplicar la extinción
def calcular_luminosidad_observada(Lbol, A_v):
    """
    Calcula la luminosidad bolométrica observada considerando la extinción Av.

    Parámetros:
    - Lbol: array o valor de luminosidades bolométricas [L_sun]
    - A_v: extinción en magnitudes

    Retorna:
    - Lbol_obs: luminosidad bolométrica observada [L_sun]
    """
    return Lbol * 10**(-0.4 * A_v)

# Inicializar listas para almacenar los resultados
time_array_L = []
L_total_Lsun = []
L_total_Lsun_obs = []  # Nueva lista para luminosidad observada

for idx, output_number in enumerate(output_numbers):
    # Calcular el tiempo correspondiente
    time_code_units = output_number * Ninter * DT
    time_days = time_code_units * unit_time / 86400  # Tiempo en días
    time_array_L.append(time_days)
    
    # Calcular \dot{M}(r, t)
    mdot_r = []
    for r0 in radii:
        mdot = compute_accretion_rate(output_number, path, r0)
        mdot_r.append(mdot)
    mdot_r = np.array(mdot_r)
    
    # Convertir \dot{M} a g/s
    mdot_r_cgs = mdot_r * (1.989e33) / (3.154e7)  # g/s
    
    # Establecer a cero los valores negativos de \dot{M}
    mdot_r_cgs = np.where(mdot_r_cgs > 0, mdot_r_cgs, 0.0)
    
    # Verificar si \dot{M} es cero en todos los radios
    if np.all(mdot_r_cgs == 0):
        print(f"At time {time_days:.1f} days, \dot{{M}} is zero or negative at all radii.")
        L_total_Lsun_value = 0.0
    else:
        # Calcular dL/dr
        dLdr = (3 * G_cgs * M_star_cgs * mdot_r_cgs) / (2 * (radii_cm ** 2))  # erg/s/cm
        
        # Integrar para obtener L_total
        L_total = np.trapz(dLdr, radii_cm)  # erg/s
        L_total_Lsun_value = L_total / L_sun_cgs  # Convertir a L_sun

    # Almacenar la luminosidad total
    L_total_Lsun.append(L_total_Lsun_value)
    
    # Calcular la luminosidad observada con extinción
    L_obs = calcular_luminosidad_observada(L_total_Lsun_value, A_V)
    L_total_Lsun_obs.append(L_obs)

# Convertir listas a arrays numpy
time_array_L = np.array(time_array_L)
L_total_Lsun_obs = np.array(L_total_Lsun_obs)

# Graficar la luminosidad bolométrica observada vs. tiempo
plt.figure(figsize=(10, 6))
plt.plot(time_array_L, L_total_Lsun_obs, marker='o', color='blue', label=f'Extinción $A_V={A_V}$ mag')
plt.xlabel('Time [days]')
plt.ylabel('Luminosidad Bolométrica Observada $L_{\\text{bol, obs}}$ [$L_\\odot$]')
plt.title(f'Luminosidad Bolométrica Total Observada $L_{{\\text{{bol, obs}}}}(t)$ con Extinción $A_V={A_V}$ mag')

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



