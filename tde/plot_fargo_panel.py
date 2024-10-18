import matplotlib.pyplot as plt
import numpy as np
import astropy 
from astropy import units as u
from astropy import constants as const
import os
import re
import glob


data_in='/Users/matias/Simulations/mi_fargo3d/outputs/tde_2d/'
path = data_in
out_folder='/Users/matias/Simulations/plotting/tde/Figures/'

snap_init=0
snap_end=100
snap_dif=1

out_name='test_'
data1_n='gasdens'
data2_n='gasenergy'
data3_n='gasvy'

gamma=1.666666666
mu=2.35 
kb_cgs=(const.k_B).cgs
R_cgs=(const.R).cgs
G_cgs=(const.G).cgs
mp_cgs=(const.m_p).cgs
R_phy_cgs=kb_cgs/(mu*mp_cgs)

print("G= {:.4e} , Kb= {:.4e}, R={:.4e}, mp={:.4e} , R_phy={:.4e}".format(G_cgs,kb_cgs,R_cgs,mp_cgs,R_phy_cgs))
# scaling units for cgs units 
unit_mass_cgs   = (const.M_sun).cgs # sun mass
unit_length_cgs = (1.0*u.AU).cgs #R0 in cm
unit_time_cgs   =  np.sqrt(unit_length_cgs**3/(G_cgs*unit_mass_cgs)) # OMEGA in s
unit_temperature= G_cgs*unit_mass_cgs*mu*mp_cgs/(kb_cgs*unit_length_cgs) #T0 in K

unit_density_cgs=unit_mass_cgs/unit_length_cgs**3 # g/cm**3
unit_surf_density_cgs = unit_mass_cgs/unit_length_cgs**2 # g/cm**2
unit_velocity_cgs = unit_length_cgs/unit_time_cgs # cm/s
unit_energy_cgs   = unit_mass_cgs*unit_length_cgs**2/unit_time_cgs**2/unit_length_cgs**3 #erg/cm**3
unit_surf_energy_cgs=unit_mass_cgs*unit_length_cgs**2/unit_time_cgs**2/unit_length_cgs**2 #erg/cm**2
unit_pressure_cgs = unit_energy_cgs

# scaling units for scale free 
unit_mass = 1.0 #
unit_length = 1.0 # 
unit_time   = 1.0 #
R_gas_mu_SF = 1.0 # gas constants R/mu

variables_par = np.genfromtxt(data_in+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Le
parametros_par, valores_par = [],[]  
for posicion in variables_par:          #
        parametros_par.append(posicion[0].decode("utf-8"))    #
        valores_par.append(posicion[1].decode("utf-8"))        
def P(parametro):
        return valores_par[parametros_par.index(parametro)] 

DT=float(P("DT"))
nintern=int(P("NINTERM"))
nx=int(P("NX"))
ny=int(P("NY"))
NX = nx
NY = ny

# Identificar el número de outputs disponibles
# Patrón para archivos gasdens*.dat con un número entero
pattern = re.compile(r'gasdens(\d+)\.dat')
# Lista de archivos que coinciden con el patrón gasdens*.dat
files = glob.glob(path + "gasdens*.dat")
# Filtrar archivos que se ajustan al patrón correcto
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
# Contar el número de archivos válidos
Ntot = len(valid_files)

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out

def cargar_energy(output_number, path):
    energy_out = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX)
    return energy_out

def cargar_temperature(output_number, path):
    dens    = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)   #no units
    energy  = np.fromfile(path + f"gasenergy{output_number}.dat").reshape(NY, NX)    #no units
    press   = (gamma-1.0)*energy 
    temp    = press/(dens*R_gas_mu_SF)
    return temp



# Calcular los valores mínimos y máximos de la densidad globalmente
vmin_gas, vmax_gas = None, None
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)*unit_density_cgs
    log_dens_out = np.log10(dens_out.value)
    if vmin_gas is None or vmax_gas is None:
        vmin_gas = log_dens_out.min()
        vmax_gas = log_dens_out.max()
    else:
        vmin_gas = min(vmin_gas, log_dens_out.min())
        vmax_gas = max(vmax_gas, log_dens_out.max())

# Calcular los valores mínimos y máximos de la densidad globalmente
vmin_temp, vmax_temp = None, None
for output_number in range(Ntot):
    temp_out = cargar_temperature(output_number, path)*unit_temperature
    log_temp_out = np.log10(temp_out.value)
    if vmin_temp is None or vmax_temp is None:
        vmin_temp = log_temp_out.min()
        vmax_temp = log_temp_out.max()
    else:
        vmin_temp = min(vmin_temp, log_temp_out.min())
        vmax_temp = max(vmax_temp, log_temp_out.max())

# Calcular los valores mínimos y máximos de la densidad globalmente
vmin_energy, vmax_energy = None, None
for output_number in range(Ntot):
    energy_out = cargar_energy(output_number, path)*unit_temperature
    log_energy_out = np.log10(energy_out.value)
    if vmin_energy is None or vmax_energy is None:
        vmin_energy = log_energy_out.min()
        vmax_energy = log_energy_out.max()
    else:
        vmin_energy = min(vmin_energy, log_energy_out.min())
        vmax_energy = max(vmax_energy, log_energy_out.max())

temp_max=vmax_temp
temp_min=vmin_temp

e_max=vmin_energy
e_min=vmin_energy

rho_max=vmax_gas
rho_min=vmin_gas

log_fixed_ranges=True

# domain files #
phi=np.loadtxt(data_in+'domain_x.dat')  # corners for phi
r=np.loadtxt(data_in+'domain_y.dat')[3:-3] # corners for r
phi_c=(phi[:-1]+phi[1:])/2.0
r_c  =(r[:-1]+r[1:])/2.0

## mesh grids ##
phi_m,r_m=np.meshgrid(phi,r)
## pas to cartesian 2d ##
x=(r_m*np.cos(phi_m)*unit_length_cgs.to(u.AU)).value
y=(r_m*np.sin(phi_m)*unit_length_cgs.to(u.AU)).value

for i in range(snap_init,snap_end+snap_dif,snap_dif):
    print("i: {}".format(i))
    # read binary data files #
    energy =np.fromfile(data_in+data2_n+str(i)+'.dat').reshape(ny,nx)
    dens=np.fromfile(data_in+data1_n+str(i)+'.dat').reshape(ny,nx)
    press=(gamma-1.0)*energy #e=Sigma*u where u=P/(sigma*(gamma-1)) is the thermal energy per unit mass, so e is the thermal energy per unit area#
    temp=press/(dens*R_gas_mu_SF)
    #temp=press/(dens)

    #print("press: {}".format(press))
    #print("energy: {}".format(energy))
    #print("unit surface energy: {}".format(unit_surf_energy_cgs))
    ### cgs ##
    energy=energy*unit_surf_energy_cgs
    dens=(dens*unit_surf_density_cgs) 
    press=press*unit_pressure_cgs  # surface pressure ???
    temp=temp*unit_temperature
    
    fig, axs=plt.subplots(nrows=1,ncols=3,figsize=(26,8))
    fig.tight_layout()
    if(log_fixed_ranges):
        pc1=axs[0].pcolormesh(x,y,np.log10(temp.value),vmin=temp_min,vmax=temp_max,cmap='inferno')
    else:
        pc1=axs[0].pcolormesh(x,y,np.log10(temp.value),vmin=temp_min,vmax=temp_max,cmap='inferno')
    axs[0].set_xlabel('x (AU)')
    axs[0].set_ylabel('y (AU)')
    axs[0].minorticks_on()
    axs[0].set_aspect('equal') 
    axs[0].title.set_text("t= {:.3e} T0 ".format(((DT*nintern*i)/(2*np.pi))))
    cb1=fig.colorbar(pc1,label='$\\log(T [K]))$')
    cb1.ax.minorticks_on()
    
    if(log_fixed_ranges):
        pc2=axs[1].pcolormesh(x,y,np.log10(energy.value),vmin=e_min,vmax=e_max)
    else:
        pc2=axs[1].pcolormesh(x,y,np.log10(energy.value),vmin=np.amin(np.log10(energy.value)),vmax=np.amax(np.log10(energy.value)))

    axs[1].set_xlabel('x (AU)')
    axs[1].set_ylabel('y (AU)')
    axs[1].minorticks_on()
    axs[1].set_aspect('equal')
    cb2=fig.colorbar(pc2,label='$\\log(e [ergs/cm^{2}] )$')
    cb2.ax.minorticks_on()
    if(log_fixed_ranges):
        pc3=axs[2].pcolormesh(x,y,np.log10(dens.value),vmin=rho_min,vmax=rho_max,cmap='plasma')
    else:
        pc3=axs[2].pcolormesh(x,y,np.log10(dens.value),vmin=np.amin(np.log10(dens.value)),vmax=np.amax(np.log10(dens.value)),cmap='plasma')

    axs[2].set_xlabel('x (AU)')
    axs[2].set_ylabel('y (AU)')
    axs[2].minorticks_on()
    axs[2].set_aspect('equal')
    cb3=fig.colorbar(pc3,label='$\\log(\\Sigma [M_{\\odot}/AU^{2}] )$')
    cb3.ax.minorticks_on()
   
    fig.tight_layout()
    #fig.colorbar(pc,label='$\\log(\\Sigma)$')
    fig.savefig(out_folder+out_name+str(i).zfill(3)+'.png')
    #plt.minorticks_on()
    #fig.savefig(data2_n+str(i)+'_FB_on.png')
    plt.close()