import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

print("\nSelecionando localização na Terra:")
loc1_lon = -35.908276*u.deg
loc1_lat =  -7.211803*u.deg
loc1 = EarthLocation.from_geodetic(lon=loc1_lon, lat=loc1_lat)
print("\nLabMet, UFCG:\n")
print(loc1)

print("\nSelecionando data para observação:")
fuso_loc1 = -3*u.hour
data_loc1 = "2022-08-12 15:00" # UTC
data_obs_loc1 = Time(data_loc1)
print("Dia e hora (UTC) - loc1:")
print(data_obs_loc1)

print("Dia e hora (loc1):")
print(data_obs_loc1 + fuso_loc1)
time_grid_loc1 = data_obs_loc1 + np.linspace(0, 365, 365)*u.d
time_grid_fuso_loc1 = (time_grid_loc1 + fuso_loc1).datetime

alt_az_loc1 = AltAz(location = loc1, obstime = time_grid_loc1)

# 1) Condições de observabilidade:
# - escolher um astro qualquer que seja do seu interesse (ex.: centro da Via Lactea, nebulosa, buraco negro, pulsar, etc)
# - Verificar se o Uirapuru é capaz de observá-lo, considerando as características de observação:
# >> az ~ 0, 40 <= el <= 80
# Calcular a elevação e hora da observação

ra1 = 17.761111*u.deg # 17 45' 40.045"
dec1 = -29.00775*u.deg
print("Astro de interesse: Sagittarius A*")
obj1 = SkyCoord(ra=ra1, dec=dec1)

print(obj1)
print("RA:")
print(obj1.ra)
print("RA em deg:")
print(obj1.ra.degree)
print("RA em horas:")
print(obj1.ra.hour)
print("Dec em string:")
print(obj1.to_string("hmsdms"))

print("\nVerificando movimento aparente:")
obj1_altaz_loc1 = obj1.transform_to(alt_az_loc1)

plt.figure(1)
plt.title("Observação Local")
plt.plot(obj1_altaz_loc1.az.degree, obj1_altaz_loc1.alt.degree, label="Sagittarius A*")
plt.plot(np.linspace(-10, 10, 100), np.linspace(40, 40, 100), label="Altura mínima")
plt.plot(np.linspace(-10, 10, 100), np.linspace(80, 80, 100), label="Altura máxima")
plt.plot(np.linspace(-50, 250, 100), np.zeros(100), label="Horizonte", c="black")
plt.ylim(-90,90)
plt.xlabel("Azimute (deg)")
plt.ylabel("Altura (deg)")
plt.legend()
plt.show()

print("\nLogo, o Sagittarius A* não é observável, pois ele nunca está na região de observação do Uirapuru.")
