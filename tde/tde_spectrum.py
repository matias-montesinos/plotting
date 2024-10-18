import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import G, c, h, k_B, M_sun, sigma_sb

# Parámetros del disco y estrella
M_star = 0.5 * M_sun  # Masa de la estrella [kg]
Mdot = 2e-6 * M_sun / u.yr  # Tasa de acreción [kg/s]
r_in = 0.1 * u.au  # Radio interno del disco
r_out = 30.0 * u.au  # Radio externo del disco
num_rings = 500  # Número de anillos en el disco
radii = np.linspace(r_in.to(u.m).value, r_out.to(u.m).value, num_rings) * u.m  # Radios en metros

# Longitudes de onda (en micrómetros)
lambda_min = 0.1 * u.um  # 0.1 micras (UV)
lambda_max = 100 * u.um  # 100 micras (IR)
wavelengths = np.logspace(np.log10(lambda_min.value), np.log10(lambda_max.value), 1000) * u.um

# Función de cuerpo negro usando astropy (B_lambda en unidades de longitud de onda)
def planck(wavelength, T):
    wavelength = wavelength.to(u.m)  # Convertimos a metros para la fórmula
    return (2 * h * c**2 / wavelength**5) / (np.exp((h * c) / (wavelength * k_B * T)) - 1)

# Temperatura en función del radio en el disco (pre-outburst)
def T_disc(r):
    r = r.to(u.m)  # Convertimos el radio a metros
    return (3 * G * M_star * Mdot / (8 * np.pi * r**3 * sigma_sb))**0.25

# Espectro del disco pre-outburst
def spectrum_pre_outburst(wavelengths):
    flux = np.zeros_like(wavelengths.value) * u.erg / (u.s * u.cm**2 * u.AA)  # Inicialización con unidades correctas
    for r in radii:
        T = T_disc(r)
        # Calcular el flujo del cuerpo negro en unidades de longitud de onda
        B_lambda = planck(wavelengths, T).to(u.erg / (u.s * u.cm**2 * u.m))  # Convertir a erg/s/cm²/m
        # Convertir de metros a Ångstrom (1 m = 10^10 Å)
        B_lambda_per_AA = B_lambda * (1e10)  # Convertimos a erg/s/cm²/Å
        # Multiplicar por el área del anillo y sumar al flujo total
        dA = 2 * np.pi * r * (radii[1] - radii[0])  # Área del anillo
        flux += B_lambda_per_AA * dA.to(u.cm**2)
    return flux

# Modelo de outburst: aumentar la temperatura en un anillo a 0.6 AU
def spectrum_outburst(wavelengths):
    flux = np.zeros_like(wavelengths.value) * u.erg / (u.s * u.cm**2 * u.AA)  # Inicialización con unidades correctas
    for r in radii:
        T = T_disc(r)
        if (0.55 * u.au < r < 0.65 * u.au):
            T *= 2.5  # Aumentar la temperatura en el anillo del outburst
        B_lambda = planck(wavelengths, T).to(u.erg / (u.s * u.cm**2 * u.m))  # Convertir a erg/s/cm²/m
        B_lambda_per_AA = B_lambda * (1e10)  # Convertir de m a Å (1 m = 10^10 Å)
        dA = 2 * np.pi * r * (radii[1] - radii[0])  # Área del anillo
        flux += B_lambda_per_AA * dA.to(u.cm**2)
    return flux

# Calcular espectros
flux_pre_outburst = spectrum_pre_outburst(wavelengths)
flux_outburst = spectrum_outburst(wavelengths)

# Crear el gráfico
plt.figure(figsize=(8, 6))
plt.plot(wavelengths.to(u.micron).value, flux_pre_outburst.value, label='Pre-outburst', color='black')
plt.plot(wavelengths.to(u.micron).value, flux_outburst.value, label='Outburst', linestyle='--', color='red')

# Configurar el gráfico
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wavelength [μm]')
plt.ylabel('Flux [erg/cm²/s/Å]')
plt.title('Ad-hoc Disc Heating Model with Outburst')
plt.legend()
plt.grid(True)

# Mostrar el gráfico
plt.show()
