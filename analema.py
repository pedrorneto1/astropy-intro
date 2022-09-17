import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun

# 2) Analema: escrever um código que faz o desenho de um analema no céu dada uma localização específica e uma hora do dia fixa.
# - fazer pelo menos 4 localizações diferentes ao meio-dia

print("\nSelecionando horário para observação...")

data_local = "2022-08-12 06:00" # NÃO ESTÁ EM UTC AINDA
data_obs = Time(data_local)


LON_KEY = "lon"
LAT_KEY = "lat"
LOC_STR_KEY = "loc"
TZ_KEY = "tz"

locs = []

locs.append({LOC_STR_KEY:"LabMet, UFCG", 
			 LON_KEY:-35.908276*u.deg,
			 LAT_KEY:-7.211803*u.deg,
			 TZ_KEY:-3*u.hour})

locs.append({LOC_STR_KEY:"Estacao Antartica Comandante Ferraz", 
			 LON_KEY:-58.391337*u.deg,
			 LAT_KEY:-62.085423*u.deg,
			 TZ_KEY:12*u.hour})
			 
locs.append({LOC_STR_KEY:"Aeroporto Internacional de Keflavik, Islandia", 
			 LON_KEY:-22.628188*u.deg,
			 LAT_KEY:63.982560*u.deg,
			 TZ_KEY:0*u.hour})

locs.append({LOC_STR_KEY:"Aeroporto Internacional do Havai", 
			 LON_KEY:-157.923837*u.deg,
			 LAT_KEY:21.334006*u.deg,
			 TZ_KEY:-10*u.hour})



az_max = 180
az_min = az_max-360

plt.figure(0)
plt.title("Analema para o Sol")

for loc in locs:

	loc_earth = EarthLocation.from_geodetic(lon=loc[LON_KEY], lat=loc[LAT_KEY])
	print("\nLocalização: {}\n".format(loc[LOC_STR_KEY]))
	print(loc_earth)
	
	# Convertendo de hora local -> UTC:
	data_utc = data_obs - loc[TZ_KEY] # Local = UTC + FH // UTC = Local - FH
	print("Dia e hora (UTC):")
	print(data_utc)
	
	time_grid = data_utc + np.linspace(0, 365, 366)*u.d
	
	print("\nBuscando Posição do Sol:")
	sun_loc = get_sun(data_utc)
	print(sun_loc)
	sun_coords = get_sun(time_grid)
	alt_az_loc = AltAz(location = loc_earth, obstime = time_grid)
	sun_altaz_loc = sun_coords.transform_to(alt_az_loc)

	az_rescale = [az if az<az_max else az-360 for az in sun_altaz_loc.az.degree]
	plt.scatter(az_rescale, sun_altaz_loc.alt.degree, label=loc[LOC_STR_KEY])
	
plt.plot(np.linspace(az_min, az_max, 100), np.zeros(100), label="Horizonte", c="black")
plt.ylim(-90,90)
plt.xlim(az_min, az_max)
plt.xlabel("Azimute (deg)")
plt.ylabel("Altura (deg)")
plt.legend()
plt.show()
