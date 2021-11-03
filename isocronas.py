# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 22:23:52 2020

@author: Juan Herrera
"""

from astropy.io.votable import parse_single_table
import numpy as np
import matplotlib.pyplot as plt
# Se lee la tabla y los valores de temperatura y luminosidad
data = parse_single_table('ic 4665-result.vot').array
temperature = data['teff_val']
luminosity = data['lum_val']
# -- En la siguiente línea, “s” es el tamaño del punto -- #
plt.scatter(temperature, luminosity, c='royalblue', s=80, \
edgecolors='k', linewidth=1.0, alpha=0.5)
# La siguiente línea lee los datos del archivo .iso
data1 = np.loadtxt('MIST_iso_5f7802c77620c.iso', comments='#').T
data2 = np.loadtxt('MIST_iso_5f7802eff1696.iso', comments='#').T
data3 = np.loadtxt('MIST_iso_5f7803b246c78.iso', comments='#').T
data4 = np.loadtxt('MIST_iso_5f78047a42361.iso', comments='#').T
data5 = np.loadtxt('MIST_iso_5f7804e6c3d38.iso', comments='#').T
LogL1 = data1[7] # Dato del logaritmo10 de Luminosidad
LogT1 = data1[10] # Dato del logaritmo10 de Temperatura
LogL2 = data2[7] # Dato del logaritmo10 de Luminosidad
LogT2 = data2[10] # Dato del logaritmo10 de Temperatura
LogL3 = data3[7] # Dato del logaritmo10 de Luminosidad
LogT3 = data3[10] # Dato del logaritmo10 de Temperatura
LogL4 = data4[7] # Dato del logaritmo10 de Luminosidad
LogT4 = data4[10] # Dato del logaritmo10 de Temperatura
LogL5 = data5[7] # Dato del logaritmo10 de Luminosidad
LogT5 = data5[10] # Dato del logaritmo10 de Temperatura
# Se grafica Luminosidad en función de la Temperatura
#plt.plot(10**LogT1, 10**LogL1, label='1 Myr', color='red')
#plt.plot(10**LogT2, 10**LogL2, label='32 Myr', color='green')
plt.plot(10**LogT3, 10**LogL3, label='100 Myr', color='yellow')
#plt.plot(10**LogT4, 10**LogL4, label='1000 Myr', color='blue')
#plt.plot(10**LogT5, 10**LogL5, label='10000 Myr', color='orange')
plt.xlabel(r'T$_{eff}$ [K]')
plt.ylabel(r'Luminosity [L$_*$/L$_\odot$]')
plt.axis([2500, 9000,0.0001, 1000000])
plt.yscale('log')
plt.gca().invert_xaxis()
plt.legend(loc='best')
plt.grid()
plt.show()