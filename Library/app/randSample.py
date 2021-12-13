import numpy as np 
import random

xmin, xmax = 0.1E6, 3.9E6
ymin, ymax = 100E5, 2500E5
zmin, zmax = 0.1E-2, 100E-2

fpout=open('rand_HPX.txt','w')
for i in range(0, int(1E6)):
    x     = (random.random())*(xmax - xmin) + xmin
    y    = (random.random())*(ymax - ymin) + ymin
    z    = (random.random())*(zmax - zmin) + zmin
    fpout.write('%f\t%f\t%f\n'%(x,y,z))

fpout.close()