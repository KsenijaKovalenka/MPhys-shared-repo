import numpy as np
import matplotlib.pyplot as plt

nk = 11

file = 'BiTeI_hr_trivial.dat'
data = np.loadtxt(file, skiprows=1, max_rows=2, dtype=int)
norb = data[0]
nr = data[1]

weights = np.loadtxt(file, skiprows=3, max_rows=int(nr/15), dtype=float).flatten()
data = np.loadtxt(file, skiprows=3 + int(nr/15))

# getting just first of array of R vectors in data 
Rvec = np.zeros((nr, 3))
for i in range(nr):
    Rvec[i] = data[18*18*i, 0:3]
# hopping coeefs
t = data[:, 5:]
# forming real hopping matrix
HR = t[:, 0] + 1j * t[:, 1]
HR = np.reshape(HR, (nr, norb, norb))

# Points defined
Gamma = np.array([0, 0, 0])
K = np.array([1/3, 1/3, 0])
M = np.array([1/2, 0, 0])
A = np.array([0, 0, 1/2])
L = np.array([1/10, 0, 1/2]) # L_x = 0.5, changed to check with smaller amout of values
H = np.array([1/3, 1/3, 1/2])

# kvec list
kList = np.array([np.linspace(-L[0], L[0], nk), 
                  np.linspace(L[1], L[1], nk), 
                  np.linspace(L[2], L[2], nk)])
kList = np.transpose(kList)

Ek = np.zeros((nk, norb))
#Hk = np.zeros((norb, norb), dtype=complex)
for ik in range(nk):
    Hk = np.zeros((norb, norb), dtype=complex)
    for ir in range(nr):
        phase = 2*np.pi*(Rvec[ir,0] * kList[ik,0] + Rvec[ir,1] * kList[ik,1] + Rvec[ir,2] * kList[ik,2])
        for iorb in range(norb):
            for jorb in range(norb):
                Hk[iorb, jorb] += HR[ir, iorb, jorb] * np.exp(-1j*phase)
    # diagonalise and put evals into corresponding energy array
    Hk = np.linalg.eig(Hk)
    Ek[ik,:] = Hk[0] 
        
for i in range(nk):
    plt.scatter(kList[:,0], Ek[:, i])