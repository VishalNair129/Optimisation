#VISHAL NAIR  NIU - 1740105
#NICO WUSTNER NIU - 1736863

#!/usr/bin/python3

# use as:
# $ python3 mapplot.py OutputExample.txt
# with matplotlib version >3.2.2 you can apply the solution from the following link:
# https://test.ocom.vn/?url=github.com/carlinmack/mplexporter/commit/739793fdda66809df80bb6b2d19c0605544f8292

import json
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import mplleaflet

# fitxer d'entrada
in_file = sys.argv[1]

# nom del fitxer de sortida: el mateix que l'entrada, canviant extensio per html
mapfile = in_file.split('.')[-2]+".html"

# convertim les dades del fitxer en un array del numpy
xy = np.genfromtxt(in_file,delimiter="|",skip_header=2)

# Dibuixem el cami amb punts vermells connectats per linies blaves
fig = plt.figure()
plt.plot(xy[:,2], xy[:,1], 'r.')
plt.plot(xy[:,2], xy[:,1], 'b')

# Creem el mapa i el guardem amb el nom guardat a mapfile
mplleaflet.show(path=mapfile)
