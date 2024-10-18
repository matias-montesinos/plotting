#Imports
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.cm     as cm

from matplotlib.patches import Ellipse
from matplotlib.colors  import ListedColormap

from matplotlib.ticker import FuncFormatter, MultipleLocator


snapshot      = 0 # int(input("Enter output: "))#   3 #100 #23
cut_x         = 128
#Crear Evolucion
output_folder = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_3d_iso/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M1_f1/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M1_f0_highres/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M1_f1_highres/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M2_f0/"

#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/juan1/"
#output_folder = "/mnt/data1/mmontesinos/paperFargo3D/M1_f1_highres/"

output_folder_0 = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_3d_isoV2/"
output_folder_1 = "/Users/matias/Simulations/mi_fargo3d/outputs/tde_3d_isoV2/"


import warnings
warnings.filterwarnings("ignore")
#----------------------------------------------------------------------#
# functions
def obtener_temp(output_folder):

    sigma  = np.fromfile(output_folder+"gasdens"+str(snapshot)+".dat").reshape(NZ,NY,NX)*unit_density       #gr/cm3
    energy = np.fromfile(output_folder+"gasenergy"+str(snapshot)+".dat").reshape(NZ,NY,NX)*unit_energy      #erg/cm3
    temp   = energy/sigma/(Rgas/mu)*(float(P("GAMMA"))-1)  # Temperatura en Kelvin
    cs     = np.sqrt(temp*float(P("GAMMA"))*Rgas/mu)/unit_length*3.154e7

    return temp

def obtener_dens(output_folder):

    sigma  = np.fromfile(output_folder+"gasdens"+str(snapshot)+".dat").reshape(NZ,NY,NX)*unit_density       #gr/cm3
    return sigma


#Grillas --------------------------------------------------------------#


def Grilla_XY():
	Phi,R  = np.meshgrid(domain_x,domain_y) #Entrega radio y co-latitud
	X,Y  = R*np.cos(Phi), R*np.sin(Phi)
	return X,Y,Phi,R
	
def Grilla_YZ():
	#Coordendas esfericas inversas : R = r sin(theta), Phi = Phi, z = r cos(theta) 
	#En este caso theta = z angulo de colatitud, mientras phi = x angulo azimutal
	
	R,T = np.meshgrid(domain_y, domain_z) #Entrega radio y azimut
	Y,Z = R*np.sin(T),R*np.cos(T)
	return Y,Z,R,T

#Cubo Radial ----------------------------------------------------------#
def radial_cube(): 
        global  y_centro
        zero_cube   = np.zeros((int(P("NZ")),int(P("NY")),int(P("NX"))))        #Crea cubo de ceros para dimecionar array
        zero_matriz = np.zeros((int(P("NY")),int(P("NX"))))                                     #Crea matriz de ceros para dimecionar array
        y_centro = 0.5*(dy[:-1] + dy[1:])                                                       #Centra para quitar 1 dimencionado
        R_centro = y_centro.reshape(int(float(P("NY"))),1)                                              #Transforma de una fila a columnas
        R_centro = (zero_matriz+R_centro)+zero_cube                                                             #Transforma columnas a matriz y de matriz a cubo
        return R_centro


#Velocida azimutal centrada en la estrella ----------------------------#
def get_fixed_vgasx(snapshot):
        R_centro = radial_cube()
        
        #Velocidad azimutal gas
        gasvx = np.fromfile(output_folder+"gasvx"+str(snapshot)+".dat").reshape(int(P("NZ")),int(P("NY")),int(P("NX"))) #Entrega: NZ matrices de NY por NX.
        
        #Velocidad azimutal planeta
        planet  = np.genfromtxt(output_folder+"planet0.dat")                            

        omega_p = planet[snapshot][-1]                                                                          #1/yr | Velocidad angular del "marco" o del planeta.
        
        #vk      = np.sqrt(G*unit_mass/((R_centro*unit_length)))                        #cm/seg # gasvt - vk < 0 esta bien porque el gas es sub-kepleriano desde la estrella
        
        gasvp = omega_p*R_centro                                                                                        #au/yr 
        gasvt = gasvx+gasvp                                                                                                     #velocidad azimutal corregida, centrada en la estrella
        
        gasvx *=1e5*unit_velocity #cm/seg | 1e5 porque unit_velocity esta en km/seg
        gasvp *=1e5*unit_velocity #cm/seg | 1e5 porque unit_velocity esta en km/seg
        gasvt *=1e5*unit_velocity #cm/seg | 1e5 porque unit_velocity esta en km/seg
        
        return gasvt

def get_altura(snapshot):
        global midplane_nz
        R_centro = radial_cube() #au
        gasvt = get_fixed_vgasx(snapshot) #au/yr
        
        density = np.fromfile(output_folder+"gasdens"+str(snapshot)+".dat").reshape(int(P("NZ")),int(P("NY")),int(P("NX"))) #Entrega: NZ matrices de NY por NX.
        energy  = np.fromfile(output_folder+"gasenergy"+str(snapshot)+".dat").reshape(int(P("NZ")),int(P("NY")),int(P("NX"))) #Entrega: NZ matrices de NY por NX.
        
        press = (float(P("GAMMA"))-1.0)*(energy*unit_energy)                                                #erg/cm3
        cs    = np.sqrt(float(P("GAMMA"))*press/(density*unit_density))                     #cm/seg
        #cs = np.sqrt(temperature*kb*gamma/(mu*mp))                                                     #cm/seg

        H = (cs/gasvt)*R_centro #au
        
        #Busca el indice NZ correspondiente al midplane:
        midplane_nz = dz.tolist()                               #convierte a lista para poder ocupar .index()
        midplane_nz = midplane_nz.index(np.pi/2.0)      #busca el valor 1.0 al dividir por pi/2 representa el misplane
        
        H = H[int(midplane_nz)-1]                                       # -1 ya que comienza a contar desde 0
        return H, cs  

def plot_altura_XY(snapshot):
	X,Y,Phi,R = Grilla_XY()
	
	H = get_altura(snapshot)

	return X,Y,H


	
def plot_altura_3D(snapshot):
	#Grilla 3D, necesita el mismo dimencionado de H, por ello x e y centro
	x_centro = 0.5*(domain_x[:-1] + domain_x[1:])						#Centra para quitar 1 dimencionado
	y_centro = 0.5*(domain_y[:-1] + domain_y[1:])						#Centra para quitar 1 dimencionado

	Phi_centro,R_centro  = np.meshgrid(x_centro,y_centro)				#Entrega radio y azimut
	X_centro  ,Y_centro  = R_centro*np.cos(Phi_centro), R_centro*np.sin(Phi_centro)
	
	H = get_altura(snapshot)
	
	return X_centro,Y_centro,H

def compute_opacity_v0(temperature_cgs,rho_cgs):
        #BELL AND LIN OPACITY 1994
        opacity = 0
        if temperature_cgs < 167.0:
                opacity = 2e-4*pow(temperature_cgs,2.0)
        else:
                if temperature_cgs < 203.0:
                        opacity = 2e16*pow(temperature_cgs,-7.0)
                else:
                        
                        if temperature_cgs < pow(2e82*rho_cgs,2./49):
                                opacity = 0.1*pow(temperature_cgs,0.5)
                        else:
                                
                                if temperature_cgs < pow(2e89*pow(rho_cgs,1./3),1./27):
                                        opacity = 2e81*pow(rho_cgs,1.0)*pow(temperature_cgs,-24.)
                                else:
                                        
                                        if (temperature_cgs < pow(1e28*pow(rho_cgs,1./3),1./7)):
                                                opacity = 1e-8*pow(rho_cgs,2./3)*pow(temperature_cgs,3.)
                                        else:
                                                
                                                if temperature_cgs < pow(1.5e56*pow(rho_cgs,2./3),0.08):
                                                        opacity = 1e-36*pow(rho_cgs,1./3)*pow(temperature_cgs,10.)
                                                else:
                                                
                                                        if temperature_cgs < pow(4.31e20*rho_cgs,2./5):
                                                                opacity = 1.5e20*pow(rho_cgs,1.)*pow(temperature_cgs,-2.5)
                                                        else:
                                                                opacity = 0.348
        return opacity 


def compute_opacity(temperature_cgs, rho_cgs):
    """
    Calcula la opacidad según Bell y Lin (1994) para arrays de temperatura y densidad.
    
    Parámetros:
    -----------
    temperature_cgs : np.ndarray
        Array de temperaturas en unidades CGS (Kelvin).
    rho_cgs : np.ndarray
        Array de densidades en unidades CGS (g/cm^3).

    Retorno:
    --------
    opacity : np.ndarray
        Array de opacidades calculadas para cada elemento de temperatura y densidad.
    """
    # Inicializamos un array de ceros con la misma forma que temperature_cgs
    opacity = np.zeros_like(temperature_cgs)

    # Condición 1: temperature < 167.0
    mask1 = temperature_cgs < 167.0
    opacity[mask1] = 2e-4 * temperature_cgs[mask1]**2.0

    # Condición 2: 167.0 <= temperature < 203.0
    mask2 = (temperature_cgs >= 167.0) & (temperature_cgs < 203.0)
    opacity[mask2] = 2e16 * temperature_cgs[mask2]**-7.0

    # Condición 3: 203.0 <= temperature < (2e82 * rho_cgs)^(2/49)
    mask3 = (temperature_cgs >= 203.0) & (temperature_cgs < np.power(2e82 * rho_cgs, 2./49))
    opacity[mask3] = 0.1 * temperature_cgs[mask3]**0.5

    # Condición 4: (2e82 * rho_cgs)^(2/49) <= temperature < (2e89 * rho_cgs^(1/3))^(1/27)
    mask4 = (temperature_cgs >= np.power(2e82 * rho_cgs, 2./49)) & (temperature_cgs < np.power(2e89 * rho_cgs**(1./3), 1./27))
    opacity[mask4] = 2e81 * rho_cgs[mask4]**1.0 * temperature_cgs[mask4]**-24.

    # Condición 5: (2e89 * rho_cgs^(1/3))^(1/27) <= temperature < (1e28 * rho_cgs^(1/3))^(1/7)
    mask5 = (temperature_cgs >= np.power(2e89 * rho_cgs**(1./3), 1./27)) & (temperature_cgs < np.power(1e28 * rho_cgs**(1./3), 1./7))
    opacity[mask5] = 1e-8 * rho_cgs[mask5]**(2./3) * temperature_cgs[mask5]**3.0

    # Condición 6: (1e28 * rho_cgs^(1/3))^(1/7) <= temperature < (1.5e56 * rho_cgs^(2/3))^0.08
    mask6 = (temperature_cgs >= np.power(1e28 * rho_cgs**(1./3), 1./7)) & (temperature_cgs < np.power(1.5e56 * rho_cgs**(2./3), 0.08))
    opacity[mask6] = 1e-36 * rho_cgs[mask6]**(1./3) * temperature_cgs[mask6]**10.

    # Condición 7: (1.5e56 * rho_cgs^(2/3))^0.08 <= temperature < (4.31e20 * rho_cgs)^(2/5)
    mask7 = (temperature_cgs >= np.power(1.5e56 * rho_cgs**(2./3), 0.08)) & (temperature_cgs < np.power(4.31e20 * rho_cgs, 2./5))
    opacity[mask7] = 1.5e20 * rho_cgs[mask7]**1.0 * temperature_cgs[mask7]**-2.5

    # Condición final (else): temperature >= (4.31e20 * rho_cgs)^(2/5)
    mask8 = temperature_cgs >= np.power(4.31e20 * rho_cgs, 2./5)
    opacity[mask8] = 0.348

    return opacity

def plot_field(field, valor):

    # Configuración de la figura
    bottom = 0.12
    left = 0.13
    top = 1. - 0.095
    right = 1. - 0.0005
    figheight = 4.88
    figwidth = 7.00

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figwidth, figheight))
    plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right)

    # Determinar los valores máximo y mínimo del campo
    print("Max field value: ", field.max())
    print("Min field value: ", field.min())

    # Ploteo del campo
    im = plt.pcolormesh(RS[:,:,cut_x] ,Z[:,:,cut_x],field,shading="gouraud",vmin=field.min(),vmax=field.max(),cmap="jet")

    # Título y etiquetas
    ax.set_title(r"Feedback off ($\kappa = " + str(valor) + r" cm^2 g^{-1}$)  Time: " + str(round(100 * unit_disktime, 2)) + " yr")
    cbar = fig.colorbar(im, ax=ax)
    #cbar.set_label(r"Temperature $K$", fontsize=14)
    cbar.set_label(r"$Log_{10}$ Effective optical depth", fontsize=14)
    ax.set_ylabel("Z [au]", fontsize=14)
    ax.set_xlabel("R [au]", fontsize=14)

    # Radios para los elipses
    radius = 0.33  # Puedes cambiarlo según condiciones
    ax.add_artist(Ellipse(xy=(1., 0), width=1 * 2 * Rhill, height=1 * 2 * Rhill, fill=False, lw=3, alpha=0.7, color="white"))  # Radio de Hill
    ax.add_artist(Ellipse(xy=(1., 0), width=radius * 2 * Rhill, height=radius * 2 * Rhill, fill=False, lw=3, alpha=0.7, color="red"))  # Radio de Ionización

    # Anotaciones
    ax.annotate('Hill radius', xy=(0.96, 0.05), xycoords='data', xytext=(0.35, 0.3), textcoords='figure fraction',
                arrowprops=dict(facecolor='white', shrink=0.1),
                horizontalalignment='right', verticalalignment='top', size=14, weight='bold', color="white")

    ax.annotate('Ionization radius', xy=(0.99, 0.02), xycoords='data', xytext=(0.3, 0.2), textcoords='figure fraction',
                arrowprops=dict(facecolor='red', shrink=0.1),
                horizontalalignment='right', verticalalignment='top', size=13, weight='bold', color="red")

    my_color = "tab:pink"
    #ax.annotate('Flow enhancement', xy=(1.03, 0.03), xycoords='data', xytext=(0.82, 0.18), textcoords='figure fraction',
    #            arrowprops=dict(facecolor=my_color, shrink=0.1),
    #            horizontalalignment='right', verticalalignment='top', size=13, weight='bold', color=my_color)

    # Mostrar la figura
    #plt.show()
    #plt.close()



#---#


#Dominios -------------------------------------------------------------#
domain_x = np.genfromtxt(output_folder+"/domain_x.dat")
domain_y = np.genfromtxt(output_folder+"/domain_y.dat")[3:-3] #[3:-3] porque dominio Y utiliza celdas "fantasmas"
domain_z = np.genfromtxt(output_folder+"/domain_z.dat")[3:-3] #[3:-3] porque dominio Z utiliza celdas "fantasmas"


print("\n")
print("Plotting model: ", output_folder)
print("Output :", snapshot)
print("\n")

#----------------------------------------------------------------------#
variables_par = np.genfromtxt(output_folder+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo variable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                     #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                          #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                              #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                 #
        
def P(parametro):
        return valores_par[parametros_par.index(parametro)]                                                                                                                                                                             #Retorna siempre str, recordad tranformar si necesita un int o float
        
#Dominios
NX = int(P("NX")); NY = int(P("NY")); NZ = int(P("NZ"))

#Unidades
G             = 6.67259e-8                                                                                              #cm3/gr/s2
mu            = 2.4                                                                                                             #gr/mol
Rgas          = 8.314472e7                                                                                              #erg/K/mol = gr*cm2/s2/K/mol
sigma_SB      = 5.67051e-5                                                                                              #erg/cm2/K4/s
c_speed       = 2.99e10                                                                                                 #cm/s
unit_length   = (1.0)*1.49598e13                                                                                #cm
unit_mass     = 1.9891e33                                                                                               #gr
unit_density  = unit_mass/unit_length**3                                                                #gr / cm3

unit_time     = np.sqrt( pow(unit_length,3.) / G / unit_mass)/3.154e7   #yr
unit_disktime = float(P("NINTERM"))*float(P("DT"))*unit_time                    #yr/snapshot      [Las unidades de tiempo del disco se toman en R=1]
unit_velocity = 1e-5*unit_length/float(unit_time*3.154e7)                               #km/s
unit_energy   = unit_mass/(unit_time*3.154e7)**2/unit_length                    #erg/cm3 = gr/s2/cm
unit_eV       = 1.60218e-12                                                                                     #erg
mproton       = 1.6726231e-24                                                                                   #gr
h_planck      = 6.6260755e-27                                                                                   #erg*s
kb            = 1.380658e-16                                                                                    #erg/K

#New Colormap
Dark_RdBu = plt.cm.RdBu_r(np.arange(plt.cm.RdBu_r.N))                                           
Dark_RdBu[:,0:3] *= 0.95                                                                                                #El color centrar no es tan blanco. 0 Negro , 1 Blanco
Dark_RdBu = ListedColormap(Dark_RdBu)                                                                   

#Dominios
print("Calculando Dominios ...")
dx = np.loadtxt(output_folder+"domain_x.dat")+np.pi/2
dy = np.loadtxt(output_folder+"domain_y.dat")[3:-3]
dz = np.loadtxt(output_folder+"domain_z.dat")[3:-3]

r = 0.5*(dy[:-1] + dy[1:]) #X-Center
p = 0.5*(dx[:-1] + dx[1:]) #X-Center
t = 0.5*(dz[:-1] + dz[1:]) #X-Center

#Grillas
print("Calcuando Grillas ...")
R,T,Phi    = np.meshgrid(r,t,p)
X,Y,Z      = R*np.cos(Phi)*np.sin(T) ,  R*np.sin(Phi)*np.sin(T) , R*np.cos(T)

RS         = np.sqrt(X**2+Y**2+Z**2)                                                                            #Radio desde la estrella.

planet  = np.genfromtxt(output_folder+"planet0.dat")                            
xp      = planet[snapshot][1] ; yp      = planet[snapshot][2]
zp      = planet[snapshot][3] ; mp      = planet[snapshot][7]
rp      = np.sqrt((xp**2)+(yp**2)+(zp**2))
Rhill   = rp*(mp/3)**(1/3)
#----------------------------------------------------------------------#
RP      = np.sqrt(((X)**2)+(Y-xp)**2+(Z)**2)
#----------------------------------------------------------------------#





temp_0 = obtener_temp(output_folder_0)
temp_1 = obtener_temp(output_folder_1)

dens_0 = obtener_dens(output_folder_0)
dens_1 = obtener_dens(output_folder_1)


#computing the optical depth
opacity_0 = 1.0
tau_0 = 0.5* opacity_0 * dens_0
tau_eff_0 = np.sqrt(3.0)/4 + 3.0*tau_0/8.0 + 1./(4.0*tau_0)

#computing the optical depth
opacity_1 = 0.01
tau_1 = 0.5* opacity_1 * dens_1
tau_eff_1 = np.sqrt(3.0)/4 + 3.0*tau_1/8.0 + 1./(4.0*tau_1)

opacity_1 = 1.0
tau_1 = 0.5* opacity_1 * dens_1

termal_0 = temp_0*dens_0
termal_1 = temp_1*dens_1

opacity = compute_opacity(temp_1, dens_1)

# Selecciona el corte adecuado en el eje X
opacity_field = opacity[:,:,cut_x]
field_0 = (termal_0[:,:,cut_x])
field_1 = (termal_1[:,:,cut_x])

field_ratio = temp_0 / temp_1
H, cs = get_altura(snapshot)
surface_dens = dens_0 / H

print(surface_dens.shape)

plot_field(np.log10(surface_dens[:,:,cut_x]), 1) #recibe un corte field[:,:,cut_x]
plot_field(np.log10(tau_eff_1[:,:,cut_x]), 0.01) #recibe un corte field[:,:,cut_x]
plt.show()
#plt.close()


print("Max field value: ", field_ratio.max())
print("Min field value: ", field_ratio.min())


# Generar el gráfico de la altura de escala del disco H en el plano medio
def plot_H(snapshot):
    # Obtener H (altura de escala) y cs (velocidad del sonido) desde la función get_altura
    H, cs = get_altura(snapshot)
    
    # Obtener la coordenada radial del plano medio
    r_midplane = radial_cube()[midplane_nz-1, :, 0]  # Selecciona el plano medio en el eje Z

    # Crear una figura y eje para el gráfico
    plt.figure(figsize=(8, 6))
    
    # Graficar H en función de la coordenada radial (r_midplane)
    plt.plot(r_midplane, H, label='Altura de escala H', color='blue')
    
    # Etiquetas y título del gráfico
    plt.xlabel('Radio (au)', fontsize=14)
    plt.ylabel('Altura de escala H (au)', fontsize=14)
    plt.title(f'Altura de escala en el plano medio para el snapshot {snapshot}', fontsize=16)
    
    # Mostrar la leyenda
    plt.legend()
    
    # Mostrar el gráfico
    plt.grid(True)
    plt.show()

# Llamar a la función para generar el gráfico
plot_H(snapshot)


def plot_H(snapshot):
    # Obtener H (altura de escala) y cs (velocidad del sonido) desde la función get_altura
    H, cs = get_altura(snapshot)
    
    # Verificar dimensiones de H
    print(f"Dimensiones de H: {H.shape}")
    
    # Obtener la coordenada radial correcta (NY) en el plano medio
    r_midplane = radial_cube()[midplane_nz-1, :, 0]  # Coordenada radial en el plano medio
    
    # Verificar dimensiones de r_midplane
    print(f"Dimensiones de r_midplane: {r_midplane.shape}")

    # Seleccionamos solamente el plano medio de H
    H_midplane = H[:, 0]  # Selecciona los valores de H en el plano medio (midplane), ajustado a la dimensión radial
    
    # Verificar que las dimensiones coincidan
    if len(r_midplane) != len(H_midplane):
        raise ValueError(f"Dimensiones incompatibles: r_midplane tiene {len(r_midplane)} elementos, pero H_midplane tiene {len(H_midplane)} elementos.")
    
    # Crear una figura y eje para el gráfico
    plt.figure(figsize=(8, 6))
    
    # Graficar H en función de la coordenada radial (r_midplane)
    plt.plot(r_midplane, H_midplane, label='Altura de escala H (Midplane)', color='blue')
    
    # Etiquetas y título del gráfico
    plt.xlabel('Radio (au)', fontsize=14)
    plt.ylabel('Altura de escala H (au)', fontsize=14)
    plt.title(f'Altura de escala en el plano medio para el snapshot {snapshot}', fontsize=16)
    
    # Mostrar la leyenda
    plt.legend()
    
    # Mostrar el gráfico
    plt.grid(True)
    plt.show()

# Llamar a la función para generar el gráfico
plot_H(snapshot)














