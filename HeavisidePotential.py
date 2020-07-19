import numpy as np
import matplotlib.pyplot as plt

N = 500
h_bar = 1
m = 0.5
L = 200 * 1e-10 / 5.29e-11
delta_x = L / N

def root(f,a):
    eps = 0.00001
    for i in range(25):
        new_a = a - f(a) * eps / (f(a+eps) - f(a))
        a = new_a
    return a

def V(i):
    x = (i - N/2)
    V0 = 0.17/13.6
    if x < - L / 2:
        return V0
    if  x > - L / 2 and x < - L / 4:
        return 2/3 * V0
    if  x > - L / 4 and x < 0:
        return 1/3 * V0
    if x >  L / 2:
        return V0
    return 0

def WaveFunc(e):
    psi = np.zeros(N)
    psi[1] = 1
    for i in range (2,N):
        psi[i] = (2*m*delta_x**2 / h_bar**2 *(v[i]-e)+2)*psi[i-1] - psi[i-2]
    return psi[N-1]

def calulation(dx = 0.0001):
    EInterval = []
    x = 0
    old = 0
    while(len(EInterval)<5):
        new = WaveFunc(x)
        if(new*old<0):
            EInterval.append(root(WaveFunc,x))
        old = new
        x += dx
    print('Eigenvalues of energy are: ',EInterval,'ev')
    return EInterval

v = list(map(V,np.arange(N)))

for e in calulation():
    psi = np.zeros(N)
    psi[1] = 0.01
    for i in range (2,N):
        psi[i] = (2*m*delta_x**2 / h_bar**2 *(v[i]-e)+2)*psi[i-1] - psi[i-2]
    plt.plot(np.arange(0,1,1/N)*200,psi*np.sqrt(1/ sum(psi * psi))+e*13.6);
plt.plot(np.arange(0,1,1/N)*200,[v[i]*13.6 for i in range(N)]);
plt.xlabel('length(Angstrom)')
plt.ylabel('potential(ev)-normal probility');
plt.show()
