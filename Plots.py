import numpy as np
import matplotlib.pyplot as plt

params = np.loadtxt(".parameters.dat")
Mb, Md, Mh = params

data = np.loadtxt("RadialVelocities.dat")
R = data[:,0]
VR = data[:,1]

bb = 0.2497
bd = 5.16
ad = 0.3105
ah = 64.3

def my_model(R, Mb, Md, Mh):
    return (np.sqrt(Mb)*R)/(R**2 + bb**2)**(3/4.) + (np.sqrt(Md)*R)/(R**2 + (bd+ad)**2)**(3/4.)  + np.sqrt(Mh)/(R**2+ah**2)**(1/4.)

VR_best = my_model(R, Mb, Md, Mh )

plt.axis([0,300,0,40])
plt.xlabel("$R$")
plt.ylabel(r"$V_R$")
plt.plot( R, VR, 'o')
plt.plot( R, VR_best, '-k', linewidth=2 )
plt.savefig("Plots.png")
