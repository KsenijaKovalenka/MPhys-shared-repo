import numpy as np
import pandas as pd
from numpy import linalg as LA
import matplotlib.pyplot as plt


nkpoints=11

klist=np.zeros((3,nkpoints))
klist[0,:]=np.linspace(-0.1,0.1,nkpoints)
klist[2,:]=0.5

hamfile="BiTeI_hr_trivial.dat"
data=np.loadtxt(hamfile,skiprows = 1,max_rows=2,dtype=int)
nband=data[0]
nrvec=data[1]
ndeg=np.loadtxt(hamfile,skiprows = 3,max_rows=int(nrvec/15),dtype=int).flatten()
data=np.loadtxt(hamfile,skiprows = int(nrvec/15)+3).reshape((nrvec,nband,nband,7))


ek=np.zeros((nband,nkpoints))
psik=np.zeros((nband,nband,nkpoints))
for ik in range(nkpoints):
	hk=np.zeros((nband,nband),dtype=complex)
	for ir in range(nrvec):
		phase=2*np.pi*data[ir,0,0,0]*klist[0,ik]+data[ir,0,0,1]*klist[1,ik]+data[ir,0,0,2]*klist[2,ik]
# =============================================================================
#       2 pi * whole phase, right???
# 		phase=2*np.pi*(data[ir,0,0,0]*klist[0,ik]+data[ir,0,0,1]*klist[1,ik]+data[ir,0,0,2]*klist[2,ik])
# =============================================================================
		for iorb in range(nband):
			for jorb in range(nband):
				hk[iorb,jorb]+=(data[ir,iorb,jorb,5]+1j*data[ir,iorb,jorb,6])*np.exp(-1j*phase)
	ek[:,ik],psik[:,:,ik]=LA.eig(hk)

for i in range(nkpoints):
    plt.scatter(klist[0,:], ek[i,:])

# =============================================================================
# plt.plot(klist[0,:], ek[11,:])
# plt.show()
# =============================================================================


#rvec=np.zeros((3,nrvec))
#hr=np.zeros((nband,nband,nrvec),dtype=complex)
#hk=np.zeros((nbnad,nband),dtype=complex)
#eigenvalus=np.zeros((nband,nkpoints)2
