import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
from PIL import Image
import glob
import re

### orbitas analiticas
# Constantes (en unidades normalizadas)
G = 1  # Constante gravitacional
M = 1  # Masa del objeto central
tolerance = 1e-10  # Tolerancia para precisión numérica

def runge_kutta_4(func, y0, t):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(1, n):
        h = t[i] - t[i-1]
        k1 = func(y[i-1], t[i-1])
        k2 = func(y[i-1] + k1 * h/2., t[i-1] + h/2.)
        k3 = func(y[i-1] + k2 * h/2., t[i-1] + h/2.)
        k4 = func(y[i-1] + k3 * h, t[i-1] + h)
        y[i] = y[i-1] + (h/6.) * (k1 + 2*k2 + 2*k3 + k4)
    return y

def derivative(state, t):
    x, y, vx, vy = state
    r = np.sqrt(x**2 + y**2)
    if r < 1e-8:  # Evitar singularidad
        r = 1e-8
    ax = -G * M * x / r**3
    ay = -G * M * y / r**3
    return np.array([vx, vy, ax, ay])

def calculate_orbit(x0, y0, v0, theta):
    vx0 = v0 * np.cos(theta)
    vy0 = v0 * np.sin(theta)
    
    r0 = np.sqrt(x0**2 + y0**2)
    
    epsilon = 0.5 * (vx0**2 + vy0**2) - G * M / r0
    h = x0 * vy0 - y0 * vx0  # Momento angular específico
    
    if abs(epsilon) < tolerance:
        # Órbita parabólica
        a = np.inf
        T_parabolica = 2 * r0 / v0
        t = np.linspace(0, 10 * T_parabolica, 10000)  # Aumentamos el tiempo y los puntos de integración
        e = 1.0
    elif epsilon < 0:
        # Órbita elíptica
        a = -G * M / (2 * epsilon)
        T = 2 * np.pi * np.sqrt(a**3 / (G * M))
        t = np.linspace(0, 2 * T, 2000)  # Dos períodos orbitales, más puntos de integración
        e = np.sqrt(1 + (2 * epsilon * h**2) / (G**2 * M**2))
    else:
        # Órbita hiperbólica
        a = -G * M / (2 * epsilon)
        v_infinito = np.sqrt(2 * epsilon)
        T_hiperbolica = r0 / v_infinito
        t = np.linspace(0, 10 * T_hiperbolica, 10000)  # Aumentamos el tiempo y los puntos de integración
        e = np.sqrt(1 + (2 * epsilon * h**2) / (G**2 * M**2))
    
    initial_state = np.array([x0, y0, vx0, vy0])
    
    solution = runge_kutta_4(derivative, initial_state, t)
    
    return solution, t, a, epsilon, e

# Parámetros iniciales
x0 = 60.  # Coordenada x inicial (en unidades astronómicas)
y0 = 10.  # Coordenada y inicial (en unidades astronómicas)
v0 = np.sqrt(2 * G * M / np.sqrt(x0**2 + y0**2))  # Velocidad inicial (velocidad de escape)
theta = np.pi  # Ángulo de lanzamiento (en radianes)

# Calcula la órbita
solution, t, a, epsilon, e = calculate_orbit(x0, y0, v0, theta)

# Extrae posiciones
x = solution[:, 0]
y = solution[:, 1]

# Gráfica de la órbita
plt.figure(figsize=(10, 10))
plt.plot(x, y)
plt.plot(0, 0, 'yo', markersize=10)
plt.title(f'Órbita del objeto (e = {e:.2f})')
plt.xlabel('Posición x (UA)')
plt.ylabel('Posición y (UA)')
plt.axis('equal')
plt.grid(True)

max_dist = 60 # max(np.max(np.abs(x)), np.max(np.abs(y)))
plt.xlim(-max_dist, max_dist)
plt.ylim(-max_dist, max_dist)
plt.plot(x0, y0, 'ro') # posición inicial del planeta
plt.arrow(x0, y0, v0 * np.cos(theta), v0 * np.sin(theta), head_width=5, head_length=5, fc='k', ec='k')  # Ajusta los parámetros para la flecha

plt.show()


# Cargar los datos
path = "/Users/matias/Simulations/mi_fargo3d/outputs/flyby2d_10MJ_NOfeel_upper/"

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

delta = 10.0

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = len(domain_x) - 1
NY = len(domain_y) - 1

dr = (Rmax - Rmin) / NY

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out

# Leer las coordenadas del planeta desde planet0.dat
def leer_coordenadas_planeta(path, snapshot, file_planet):
    planet_file0 = path + file_planet
    planet0 = np.genfromtxt(planet_file0)
    xp = planet0[snapshot][1]  # Coordenada x del planeta
    yp = planet0[snapshot][2]  # Coordenada y del planeta

    return xp, yp #, xp1, yp1

def Grilla_XY():
    R = 0.5 * (domain_y[1:] + domain_y[:-1])
    Phi = 0.5 * (domain_x[1:] + domain_x[:-1])
    P, R = np.meshgrid(Phi, R)
    X, Y = R * np.cos(P), R * np.sin(P)
    return X, Y

X, Y = Grilla_XY()

# Establecer límites fijos para los ejes
x_limits = (-(Rmax+delta), (Rmax+delta))
y_limits = (-(Rmax+delta), (Rmax+delta))

# Crear la carpeta para guardar los PNGs si no existe
output_dir = os.path.join(path, "gas_png")
os.makedirs(output_dir, exist_ok=True)


# Calcular vmin y vmax globales
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


### graficos de tests:
# Solicitar al usuario el snapshot a usar
snapshot = int(input("Ingrese el número de snapshot: "))
dens_out = cargar_densidad(snapshot, path)
xp0, yp0 = leer_coordenadas_planeta(path, snapshot, "planet0.dat")
plt.figure(figsize=(8, 8))
mesh = plt.pcolormesh(X, Y, np.log10(dens_out), cmap='jet', shading='nearest', vmin = vmin, vmax=vmax)
plt.scatter(xp0, yp0, c='teal', s=50, edgecolor='black', marker='x', label='Planet Position 1')
cbar = plt.colorbar(mesh, label='Log Gas Density [gr / cm$^2$]')
#cbar.set_clim(-4.72, -3.24)  # Ajusta los límites de la barra de colores
plt.xlabel('X [AU]', fontsize=14)
plt.ylabel('Y [AU]', fontsize=14)
plt.title('Gas Density Distribution', fontsize=16)
plt.xlim(x_limits)
plt.ylim(y_limits)
plt.gca().set_aspect('equal', adjustable='datalim')
#plt.axis('equal')
plt.legend()
plt.tight_layout()
# grafico de la orbita
plt.plot(x, y)
plt.plot(0, 0, 'yo', markersize=10)
plt.plot(x0, y0, 'ro') # posición inicial del planeta
plt.arrow(x0, y0, v0 * np.cos(theta), v0 * np.sin(theta), head_width=5, head_length=5, fc='k', ec='k')  # Ajusta los parámetros para la flecha
plt.show()



# Cargar el último archivo de salida para calcular vmin y vmax
ultimo_output = Ntot - 1
dens_out = cargar_densidad(ultimo_output, path)
log_dens_out = np.log10(dens_out)

#vmin = log_dens_out.min()
#vmax = log_dens_out.max()

# Generar los plots para cada output
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)
    xp0, yp0 = leer_coordenadas_planeta(path, output_number, "planet0.dat")
    plt.figure(figsize=(8, 8))
    mesh = plt.pcolormesh(X, Y, np.log10(dens_out), cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
    plt.scatter(xp0, yp0, c='teal', s=50, edgecolor='black', marker='x', label='Planet Position 1')
    plt.colorbar(mesh, label='Log Gas Density [gr / cm$^2$]')
    plt.xlabel('X [AU]', fontsize=14)
    plt.ylabel('Y [AU]', fontsize=14)
    plt.title(f'Gas Density Distribution - Output {output_number}', fontsize=16)
    plt.xlim(x_limits)
    plt.ylim(y_limits)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.legend()
    plt.tight_layout()

    # grafico de la orbita
    plt.plot(x, y)
    plt.plot(0, 0, 'yo', markersize=10)
    plt.plot(x0, y0, 'ro') # posición inicial del planeta
    plt.arrow(x0, y0, v0 * np.cos(theta), v0 * np.sin(theta), head_width=5, head_length=5, fc='k', ec='k')  # Ajusta los parámetros para la flecha
    # Guardar el PNG en la carpeta especificada
    file_path = os.path.join(output_dir, f"output_{output_number:04d}.png")
    plt.savefig(file_path)
    plt.close()

### Animar los PNGs
# Ruta donde están almacenados los archivos PNG
png_dir = os.path.join(path, "gas_png")

# Obtener una lista de todos los archivos PNG
png_files = [os.path.join(png_dir, f"output_{i:04d}.png") for i in range(Ntot)]

# Crear la figura para la animación
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

# Ocultar los ejes y eliminar los márgenes
ax.axis('off')
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

# Función para cargar las imágenes
def load_image(file):
    img = Image.open(file)
    return [plt.imshow(img, animated=True)]

# Crear la animación
frames = [load_image(file) for file in png_files]
ani = animation.ArtistAnimation(fig, frames, interval=100, repeat_delay=1000, blit=True)

# Guardar la animación en un archivo MP4
output_file = os.path.join(png_dir, 'gas_density_animation_from_png.mp4')
ani.save(output_file, fps=5, extra_args=['-vcodec', 'libx264'])

plt.show()
