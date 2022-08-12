import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun

# 2) Analema: escrever um código que faz o desenho de um analema no céu dada uma localização específica e uma hora do dia fixa.
# - fazer pelo menos 4 localizações diferentes ao meio-dia

print("\nSelecionando localizações na Terra:")

loc1_lon = -35.908276*u.deg
loc1_lat =  -7.211803*u.deg
loc1 = EarthLocation.from_geodetic(lon=loc1_lon, lat=loc1_lat)
print("\nLabMet, UFCG:\n")
print(loc1)

loc2_lon = -58.391337*u.deg
loc2_lat =  -62.085423*u.deg
loc2 = EarthLocation.from_geodetic(lon=loc2_lon, lat=loc2_lat)
print("\nEstacao Antartica Comandante Ferraz:\n")
print(loc2)

loc3_lon = -22.628188*u.deg
loc3_lat =  63.982560*u.deg
loc3 = EarthLocation.from_geodetic(lon=loc3_lon, lat=loc3_lat)
print("\nAeroporto Internacional de Keflavik, Islandia:\n")
print(loc3)

loc4_lon = -157.923837*u.deg
loc4_lat =  21.334006*u.deg
loc4 = EarthLocation.from_geodetic(lon=loc4_lon, lat=loc4_lat)
print("\nAeroporto Internacional do Havai:\n")
print(loc4)

print("\nSelecionando datas para observação:")

fuso_loc1 = -3*u.hour
data_loc1 = "2022-08-12 15:00" # UTC
data_obs_loc1 = Time(data_loc1)
print("Dia e hora (UTC) - loc1:")
print(data_obs_loc1)

fuso_loc2 = 12*u.hour
data_loc2 = "2022-08-12 00:00" # UTC
data_obs_loc2 = Time(data_loc2)
print("Dia e hora (UTC) - loc2:")
print(data_obs_loc2)

fuso_loc3 = 0*u.hour
data_loc3 = "2022-08-12 12:00" # UTC
data_obs_loc3 = Time(data_loc3)
print("Dia e hora (UTC) - loc3:")
print(data_obs_loc3)

fuso_loc4 = -10*u.hour
data_loc4 = "2022-08-12 22:00" # UTC
data_obs_loc4 = Time(data_loc4)
print("Dia e hora (UTC) - loc4:")
print(data_obs_loc4)

# loc1
print("Dia e hora (loc1):")
print(data_obs_loc1 + fuso_loc1)
time_grid_loc1 = data_obs_loc1 + np.linspace(0, 365, 365)*u.d
time_grid_fuso_loc1 = (time_grid_loc1 + fuso_loc1).datetime

# loc2
print("Dia e hora (loc2):")
print(data_obs_loc2 + fuso_loc2)
time_grid_loc2 = data_obs_loc2 + np.linspace(0, 365, 365)*u.d
time_grid_fuso_loc2 = (time_grid_loc2 + fuso_loc2).datetime

# loc3
print("Dia e hora (loc3):")
print(data_obs_loc3 + fuso_loc3)
time_grid_loc3 = data_obs_loc3 + np.linspace(0, 365, 365)*u.d
time_grid_fuso_loc3 = (time_grid_loc3 + fuso_loc3).datetime

# loc4
print("Dia e hora (loc4):")
print(data_obs_loc4 + fuso_loc4)
time_grid_loc4 = data_obs_loc4 + np.linspace(0, 365, 365)*u.d
time_grid_fuso_loc4 = (time_grid_loc4 + fuso_loc4).datetime

print("\nBuscando Posição do Sol:")

# loc1
sun_loc1 = get_sun(data_obs_loc1)
print(sun_loc1)
sun_coords = get_sun(time_grid_loc1)
alt_az_loc1 = AltAz(location = loc1, obstime = time_grid_loc1)
sun_altaz_loc1 = sun_coords.transform_to(alt_az_loc1)

# loc2
sun_loc2 = get_sun(data_obs_loc2)
print(sun_loc2)
sun_coords = get_sun(time_grid_loc2)
alt_az_loc2 = AltAz(location = loc2, obstime = time_grid_loc2)
sun_altaz_loc2 = sun_coords.transform_to(alt_az_loc2)

# loc3
sun_loc3 = get_sun(data_obs_loc3)
print(sun_loc3)
sun_coords = get_sun(time_grid_loc3)
alt_az_loc3 = AltAz(location = loc3, obstime = time_grid_loc3)
sun_altaz_loc3 = sun_coords.transform_to(alt_az_loc3)

# loc4
sun_loc4 = get_sun(data_obs_loc4)
print(sun_loc4)
sun_coords = get_sun(time_grid_loc4)
alt_az_loc4 = AltAz(location = loc4, obstime = time_grid_loc4)
sun_altaz_loc4 = sun_coords.transform_to(alt_az_loc4)

plt.figure(0)
plt.title("Observação Local para o Sol - loc1")
az_max = 180
az_min = az_max-360
az_rescale_loc1 = [az if az<az_max else az-360 for az in sun_altaz_loc1.az.degree]
az_rescale_loc2 = [az if az<az_max else az-360 for az in sun_altaz_loc2.az.degree]
az_rescale_loc3 = [az if az<az_max else az-360 for az in sun_altaz_loc3.az.degree]
az_rescale_loc4 = [az if az<az_max else az-360 for az in sun_altaz_loc4.az.degree]

plt.plot(az_rescale_loc1, sun_altaz_loc1.alt.degree, label="Sol - LabMet", c="yellow")
plt.plot(az_rescale_loc2, sun_altaz_loc2.alt.degree, label="Sol - Antartica", c="orange")
plt.plot(az_rescale_loc3, sun_altaz_loc3.alt.degree, label="Sol - Islandia", c="red")
plt.plot(az_rescale_loc4, sun_altaz_loc4.alt.degree, label="Sol - Havai", c="brown")
plt.plot(np.linspace(az_min, az_max, 100), np.zeros(100), label="Horizonte", c="black")
plt.ylim(-90,90)
plt.xlim(az_min, az_max)
plt.xlabel("Azimute (deg)")
plt.ylabel("Altura (deg)")
plt.legend()
plt.show()
