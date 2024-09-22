#Imports
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.cm     as cm

from matplotlib.patches import Ellipse
from matplotlib.colors  import ListedColormap

from matplotlib.ticker import FuncFormatter, MultipleLocator


import warnings
warnings.filterwarnings("ignore")
#-----------------------------------------------------------------------#
def calculate_X_final(temp):
    """

    """

    # Calcular exponente y exponencial
    exponente = np.float128(-(13.6 * unit_eV) / (kb * temp))
    exponencial = np.float128(np.exp(exponente))

    # Primer término
    termino_gamma = np.float128((2 * np.pi * mproton * kb * temp / h_planck**2)**(3/2))
    Gamma = np.float128(termino_gamma * exponencial)

    # Calcular X_final
    n_numero = np.float128(1.0)
    negative = np.float128(-1.0)

    X_temrino = np.float128(negative * Gamma / n_numero)
    X_square_1 = np.float128(Gamma**2 / n_numero**2)
    X_square_2 = np.float128(4.0 * Gamma / n_numero)

    X_square = np.float128(np.sqrt(X_square_1 + X_square_2))
    X_final = (X_temrino + X_square) / 2.0

    return X_final


def plot_X_final(X_final, cut_z=31, cut_y=67, cut_x=128, show_cortes=False, log_scale=False, cmap="cividis", save_fig=False, snapshot=0):
    """
    Plotea X_final en tres planos: XY, YZ y XZ.

    Parámetros:
    -----------
    X_final : np.ndarray
        Campo 3D de fracción de ionización.
    cut_z, cut_y, cut_x : int, opcional
        Cortes en los planos Z, Y y X respectivamente. Por defecto son (31, 67, 128).
    show_cortes : bool, opcional
        Muestra los cortes en la grilla. Por defecto es False.
    log_scale : bool, opcional
        Aplica escala logarítmica al colormap. Por defecto es False.
    cmap : str, opcional
        Colormap a utilizar. Por defecto es "viridis".
    save_fig : bool, opcional
        Si True, guarda la figura en lugar de mostrarla. Por defecto es False.
    snapshot : int, opcional
        Número del snapshot, utilizado para guardar la figura. Por defecto es 0.
    """

    # Escala logarítmica
    if log_scale:
        X_final = np.log10(X_final)

    # Definir los límites de color específicos para cada corte
    vmin_xy_back_2, vmax_xy_back_2 = -20, 0
    vmin_yz_back_2, vmax_yz_back_2 = -20, 0
    vmin_xz_back_2, vmax_xz_back_2 = -20, 0

    # Opciones de la figura
    bottom = 0.13
    left   = 0.05
    top    = 1.0 - 0.12
    right  = 1.0 - 0.07
    wspace = 0.4
    hspace = 0.0

    figheight = 4
    figwidth  = 16.8

    # Crear figura y ejes
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(figwidth, figheight))
    plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=wspace, hspace=hspace)
    ax = axes.flatten()

    # Alpha para transparencia si show_cortes es True
    alpha = 0.6 if show_cortes else 1.0

    # Loop para generar las tres vistas: XY, YZ, XZ
    for k in range(0, 3):
        if k == 0:  # Plano XY
            if show_cortes:
                ax[k].scatter(X[cut_z, :, cut_x], Y[cut_z, :, cut_x], s=2, color="black", marker='D')
                ax[k].scatter(X[cut_z, cut_y, :], Y[cut_z, cut_y, :], s=2, color="black", marker='D')

            ax[k].text(0.1, 0.9, "XY", transform=ax[k].transAxes, bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.4, 'pad': 0.2}, fontsize=14)

            # Cambia los límites de color aquí
            im = ax[k].pcolormesh(X[cut_z, :, :], Y[cut_z, :, :], X_final[cut_z, :, :], shading="gouraud", vmin=vmin_xy_back_2, vmax=vmax_xy_back_2, cmap=cmap, alpha=alpha)

            # Dibujar el Radio de Hill y R_eff_X
            ax[k].add_artist(Ellipse(xy=(yp, xp), width=2*Rhill, height=2*Rhill, fill=False, lw=2, alpha=0.7, color="white"))
            ax[k].add_artist(Ellipse(xy=(yp, xp), width=2*R_eff_X, height=2*R_eff_X, fill=False, lw=2, alpha=0.7, color="red"))

            ax[k].set_ylabel("au", fontsize=15)

            # Agregar la barra de colores
            cbar = fig.colorbar(im, ax=ax[k], orientation="vertical")
            cbar.set_label(r"X_final", fontsize=12)

        elif k == 1:  # Plano YZ
            if show_cortes:
                ax[k].scatter(RS[cut_z, :, cut_x], Z[cut_z, :, cut_x], s=2, color="black", marker='D')
                ax[k].scatter(RS[:, cut_y, cut_x], Z[:, cut_y, cut_x], s=2, color="black", marker='D')

            ax[k].text(0.1, 0.9, "YZ", transform=ax[k].transAxes, bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.4, 'pad': 0.2}, fontsize=14)

            # Cambia los límites de color aquí
            im = ax[k].pcolormesh(RS[:, :, cut_x], Z[:, :, cut_x], X_final[:, :, cut_x], shading="gouraud", vmin=vmin_yz_back_2, vmax=vmax_yz_back_2, cmap=cmap, alpha=alpha)

            # Dibujar el Radio de Hill y R_eff_X
            ax[k].add_artist(Ellipse(xy=(xp, zp), width=2*Rhill, height=2*Rhill, fill=False, lw=2, alpha=0.7, color="white"))
            ax[k].add_artist(Ellipse(xy=(xp, zp), width=2*R_eff_X, height=2*R_eff_X, fill=False, lw=2, alpha=0.7, color="red"))

            ax[k].set_title("Ionization region")

            # Agregar la barra de colores
            cbar = fig.colorbar(im, ax=ax[k], orientation="vertical")
            cbar.set_label(r"X_final", fontsize=12)

        elif k == 2:  # Plano XZ
            if show_cortes:
                ax[k].scatter(X[cut_z, cut_y, :], Z[cut_z, cut_y, :], s=2, color="black", marker='D')
                ax[k].scatter(X[:, cut_y, cut_x], Z[:, cut_y, cut_x], s=2, color="black", marker='D')

            ax[k].text(0.1, 0.9, "XZ", transform=ax[k].transAxes, bbox={'boxstyle': 'round', 'facecolor': 'white', 'alpha': 0.4, 'pad': 0.2}, fontsize=14)

            # Cambia los límites de color aquí
            im = ax[k].pcolormesh(X[:, cut_y, :], Z[:, cut_y, :], X_final[:, cut_y, :], shading="gouraud", vmin=vmin_xz_back_2, vmax=vmax_xz_back_2, cmap=cmap, alpha=alpha)

            # Dibujar el Radio de Hill y R_eff_X
            ax[k].add_artist(Ellipse(xy=(yp, zp), width=2*Rhill, height=2*Rhill, fill=False, lw=2, alpha=0.7, color="white"))
            ax[k].add_artist(Ellipse(xy=(yp, zp), width=2*R_eff_X, height=2*R_eff_X, fill=False, lw=2, alpha=0.7, color="red"))

            # Ajustar formato de pi para el eje azimutal
            ax[k].xaxis.set_major_formatter(FuncFormatter(lambda val, pos: '{:.0g}$\pi$'.format(val/np.pi) if val != 0 else '0'))
            ax[k].xaxis.set_major_locator(MultipleLocator(base=np.pi))

            # Agregar la barra de colores
            cbar = fig.colorbar(im, ax=ax[k], orientation="vertical")
            cbar.set_label(r"X_final", fontsize=12)

        # Configuración final del gráfico
        ax[k].tick_params(axis="both", which="major", direction="out", length=3, width=1, color="black")
        ax[k].set_xlabel("au", fontsize=15)

    # Guardar o mostrar figura
    if save_fig:
        plt.savefig(f"X_final_snapshot_{snapshot}.png", dpi=300)
    else:
        plt.show()

    plt.close()


def get_r_eff(mapa,mapa_radios,condition):
        r_effs = []
        pos_eff = []
        for i in range(0,NZ):
                for j in range(0,NY):
                        for k in range(0,NX):
                                if mapa[i][j][k] >= condition:
                                        r_effs.append(mapa_radios[i][j][k])
                                        pos_eff.append([i,j,k])
                                        #print("r: ",mapa_radios[i][j][k]," au","X: ",mapa[i][j][k],"indices" ,i, j, k)
        
        idx = r_effs.index(np.max(r_effs))
        print(idx, pos_eff[idx], r_effs[idx])
        #Toma el maximo valor de r effs
        return np.max(r_effs)


#Output
#output_folder = "/mnt/data1/mmontesinos/2021/fargo3d/model_f1/"
#output_folder = "/mnt/data1/mmontesinos/2021/fargo3d/paper_f0/"
#output_folder = "/mnt/data1/mmontesinos/2022/f3/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M1_f1/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M2_f0/"

output_folder_0 =  "/mnt/data1/mmontesinos/paperFargo3D/M1_f1/"
output_folder_1 =  "/mnt/data1/mmontesinos/paperFargo3D/M1_f1_k2/"

print(output_folder_1)

snapshot      = int(input("Enter output: ")) 

#Parametros
variables_par = np.genfromtxt(output_folder_0+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo variable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                     #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                          #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                              #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                 #
        
def P(parametro):
        return valores_par[parametros_par.index(parametro)]                                                                                                                                                                             #Retorna siempre str, recordad tranformar si necesita un int o float
        
#Dominios
NX = int(P("NX")); NY = int(P("NY")); NZ = int(P("NZ"))
print(NX, NY, NZ)

#Unidades
R_scale = 1.0
G             = 6.67259e-8                                                                                              #cm3/gr/s2
mu            = 2.4                                                                                                             #gr/mol
Rgas          = 8.314472e7                                                                                              #erg/K/mol = gr*cm2/s2/K/mol
sigma_SB      = 5.67051e-5                                                                                              #erg/cm2/K4/s
mproton       = 1.6726231e-24                                                                                   #gr
h_planck      = 6.6260755e-27                                                                                   #erg*s
kb            = 1.380658e-16                                                                                    #erg/K
unit_eV       = 1.60218e-12                                                                                     #erg
unit_length   = (1.0)*1.49598e13                                                                                #       cm
unit_mass     = 1.9891e33                                                                                               #gr

unit_density  = unit_mass/(R_scale*unit_length)**3                                                              #gr / cm3
unit_time     = np.sqrt( pow((R_scale*unit_length),3.) / G / unit_mass)/3.154e7 #yr
unit_disktime = float(P("NINTERM"))*float(P("DT"))*unit_time                    #yr/snapshot      [Las unidades de tiempo del disco se toman en R=1]
unit_velocity = 1e-5*(unit_length*R_scale)/float(unit_time*3.154e7)                             #km/s
unit_energy   = unit_mass/(unit_time*3.154e7)**2/(R_scale*unit_length)                  #erg/cm3 = gr/s2/cm


#Planet
print("Leyendo Planeta ...")
#planet  = np.genfromtxt(output_folder+"planet0.dat")                           
planet  = np.genfromtxt(output_folder_0+"bigplanet0.dat") 
xp      = planet[snapshot][1] ; yp      = planet[snapshot][2]*R_scale
zp      = planet[snapshot][3] ; mp      = planet[snapshot][7]
rp      = np.sqrt((xp**2)+(yp**2)+(zp**2))

Rhill   = rp*(mp/3)**(1./3)

#Dominios
print("Calculando Dominios ...")
dx = np.loadtxt(output_folder_0+"domain_x.dat")+np.pi/2
dy = np.loadtxt(output_folder_0+"domain_y.dat")[3:-3]*R_scale
dz = np.loadtxt(output_folder_0+"domain_z.dat")[3:-3]

r = 0.5*(dy[:-1] + dy[1:]) #X-Center
p = 0.5*(dx[:-1] + dx[1:]) #X-Center
t = 0.5*(dz[:-1] + dz[1:]) #X-Center

#Grillas
print("Calcuando Grillas ...")
R,T,Phi = np.meshgrid(r,t,p)
X,Y,Z = R*np.cos(Phi)*np.sin(T) ,  R*np.sin(Phi)*np.sin(T) , R*np.cos(T)
RS    = np.sqrt(X**2+Y**2+Z**2)   
RP    = np.sqrt(((X)**2)+(Y-xp)**2+(Z)**2)                                                                           #Radio desde la estrella.

#Mapas
print("Leyendo Densidad, Temperaturam, Energia, Cs ...")
#sigma  = np.fromfile(output_folder_0+"gasdens"+str(snapshot)+".dat").reshape(NZ,NY,NX)*unit_density       #gr/cm3
#energy = np.fromfile(output_folder_0+"gasenergy"+str(snapshot)+".dat").reshape(NZ,NY,NX)*unit_energy      #erg/cm3
#temp   = energy/sigma/(Rgas/mu)*(float(P("GAMMA"))-1)  

def obtener_campos(output_folder, snapshot):

    # Cargar densidad y energía
    sigma = np.fromfile(output_folder + "gasdens" + str(snapshot) + ".dat").reshape(NZ, NY, NX) * unit_density  # gr/cm3
    energy = np.fromfile(output_folder + "gasenergy" + str(snapshot) + ".dat").reshape(NZ, NY, NX) * unit_energy  # erg/cm3
    
    # Calcular el campo de temperatura
    temp = energy / sigma / (Rgas / mu) * (float(P("GAMMA")) - 1)  # Temperatura en K
    
    return temp

temp_0 = obtener_campos(output_folder_0, snapshot)
temp_1 = obtener_campos(output_folder_1, snapshot)

X_final_0 = calculate_X_final(temp_0)
X_final_1 = calculate_X_final(temp_1)

X_ratio = X_final_0 / X_final_1

# Mostrar el valor máximo de X_final
print(X_final_0.max())


X_condition = 0.5 #1e-6
R_eff_X     =  get_r_eff(X_final_1,RP,X_condition)
print(X_final_0/X_final_1)


# Llamada a la función para plotear X_final
print("R_eff_X:",round(R_eff_X,5)," au","for X>=",X_condition)
print("\n")
print("Rhill:", Rhill)
print("R_eff_X / Rhill:", round(R_eff_X,5)/Rhill)

# Si quieres aplicar escala logarítmica:
plot_X_final(X_final_0, show_cortes=False, log_scale=True, cmap="viridis", save_fig=False, snapshot=snapshot)

# O sin escala logarítmica:
#plot_X_final(X_final, show_cortes=False, log_scale=False, cmap="plasma", save_fig=False, snapshot=snapshot)

