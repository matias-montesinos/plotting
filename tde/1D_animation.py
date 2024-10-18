import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Ruta donde está almacenado el archivo planet0.dat
path = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d_ad/"

variables_par = np.genfromtxt(path + "/variables.par", dtype={'names': ("parametros", "valores"), 'formats': ("|S30", "|S300")}).tolist()
parametros_par, valores_par = [], []
for posicion in variables_par:
    parametros_par.append(posicion[0].decode("utf-8"))
    valores_par.append(posicion[1].decode("utf-8"))

def P(parametro):
    return valores_par[parametros_par.index(parametro)] 

# UNITS
gamma = 1.6666
mu = 2.35
kb = 1.380650424e-16  # ergs
R_gas = 8.314472e7  # erg/K/mol 
G = 6.67259e-8  # cm**3 g**-1 s**-2
mp = 1.672623099e-24  # g
unit_mass = 1.9891e33  # g
unit_length = 1.49598e13  # cm

unit_density = (unit_mass) / (unit_length**3)  # g/cm**3
unit_surf_density = (unit_mass) / (unit_length**2)  # g/cm**2
unit_time = np.sqrt(pow(unit_length, 3.) / G / unit_mass) / 3.154e7  # yr 

unit_energy = unit_mass / (unit_time * 3.154e7)**2 / (unit_length)  # erg/cm3 = gr/s2/cm
unit_temperature = ((G * mp * mu) / (kb)) * (unit_mass / (unit_length))  # K

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out * unit_surf_density

# Función para cargar la temperatura de gas
def cargar_temperature(output_number, path):
    dens = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX) * unit_density
    energy = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX) * unit_energy
    press = (gamma - 1.0) * energy
    temp = mu * press / (dens * R_gas)
    return temp

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]  # Mantener los límites fantasma

NX = int(P("NX"))  # len(domain_x) - 1
NY = int(P("NY"))  # len(domain_y) - 1
Ninter = int(P("NINTERM"))
DT = float(P("DT"))

# Calcular el número total de snapshots
Ntot = 35  # Ajustar según el número real de snapshots disponibles

# Configuración de la figura
fig, ax1 = plt.subplots(figsize=(10, 6))

# Crear un segundo eje y para la temperatura
ax2 = ax1.twinx()

# Limitar los ejes verticales de densidad y temperatura (se pueden ajustar los valores)
# Inicializar variables para los límites
vmin_density, vmax_density = np.inf, -np.inf
vmin_temp, vmax_temp = np.inf, -np.inf

# Calcular los límites de los ejes para toda la animación
for snapshot in range(Ntot):
    densidad = cargar_densidad(snapshot, path)
    temperatura = cargar_temperature(snapshot, path)
    densidad_promedio = np.mean(densidad, axis=1)
    temperatura_promedio = np.mean(temperatura, axis=1)

    vmin_density = min(vmin_density, densidad_promedio.min())
    vmax_density = max(vmax_density, densidad_promedio.max())

    vmin_temp = min(vmin_temp, temperatura_promedio.min())
    vmax_temp = max(vmax_temp, temperatura_promedio.max())

# Añadir márgenes para que los gráficos no estén ajustados
margen_density = 0.1 * (vmax_density - vmin_density)
margen_temp = 0.1 * (vmax_temp - vmin_temp)

vmin_density -= margen_density
vmax_density += margen_density

vmin_temp -= margen_temp
vmax_temp += margen_temp

dens_min, dens_max = vmin_density, vmax_density
temp_min, temp_max = vmin_temp, vmax_temp


# Inicializar las líneas que se actualizarán
line1, = ax1.plot([], [], color='#2b6ca3', linestyle='-', lw=3, alpha=0.8, label="Density")
line2, = ax2.plot([], [], color='#d35400', linestyle='-', lw=3, alpha=0.8, label="Temperature")

# Función de inicialización para la animación
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2

# Función de actualización de cada frame
def actualizar(frame):
    global dens_min, dens_max, temp_min, temp_max
    # Tiempo de evolución
    tiempo_evolucion = Ninter * DT * frame

    # Cargar los datos de densidad y temperatura del frame actual
    densidad = cargar_densidad(frame, path)
    temperatura = cargar_temperature(frame, path)

    # Tomar el promedio azimutal para obtener perfiles radiales
    densidad_promedio = np.mean(densidad, axis=1)
    temperatura_promedio = np.mean(temperatura, axis=1)

    # Asegurarse de que la grilla y los datos tengan la misma longitud
    domain_y_cut = domain_y[:NY]

    # Actualizar los datos de las líneas
    line1.set_data(domain_y_cut, densidad_promedio)
    line2.set_data(domain_y_cut, temperatura_promedio)

    # Mantener los límites fijos para los ejes
    ax1.set_xlim(domain_y_cut.min(), domain_y_cut.max())

    # Variable para activar o desactivar el uso de la escala logarítmica
    use_log_scale = True  # Cambia a False para usar escala lineal

    if use_log_scale:
        # Usar escala logarítmica
        dens_min = max(dens_min, 5e1)  # Evitar valores muy pequeños
        temp_min = max(temp_min, 1e1)

        ax1.set_yscale('log')
        ax2.set_yscale('log')
    else:
        # Usar escala lineal
        ax1.set_yscale('linear')
        ax2.set_yscale('linear')

    # Mantener los límites
    ax1.set_ylim(dens_min, dens_max)
    ax2.set_ylim(temp_min, temp_max)

    # Actualizar etiquetas y título
    ax1.set_xlabel('Radius [AU]', fontsize=16)
    ax1.set_ylabel('Density [$g/cm^2$]', color='#2b6ca3', fontsize=16)
    ax2.set_ylabel('Temperature [K]', color='#d35400', fontsize=16)

    ax1.tick_params(axis='both', labelsize=14)
    ax1.tick_params(axis='y', labelcolor='#2b6ca3')
    ax2.tick_params(axis='y', labelsize=14)
    ax2.tick_params(axis='y', labelcolor='#d35400')

    # Actualizar el título con el frame y tiempo de evolución
    fig.suptitle(f'Frame: {frame}, Evolution time: {tiempo_evolucion:.2e} yr', fontsize=16)

    return line1, line2


# Crear la animación
ani = animation.FuncAnimation(fig, actualizar, frames=Ntot, init_func=init, interval=200, repeat=False)

# Guardar la animación en un archivo MP4
output_file = os.path.join(path, '1D_dens_temp_animation.mp4')
ani.save(output_file, writer='ffmpeg', fps=10)

# Mostrar la animación en pantalla
plt.show()
