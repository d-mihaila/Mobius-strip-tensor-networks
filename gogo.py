from tkinter import Y
import numpy as np
from scipy.sparse import diags
from pytenet import OpChain
from pytenet import MPO
from pytenet import minimization
from pytenet import MPS
import matplotlib.pyplot as plt


# constants 
a = np.array([[0,1],[0,0]])
aH = a.conjugate().T
I = np.identity(2)

c = 1
ECb = 0.2*c
ECj = 100*c
EJ = 19*c

n0 = (ECj/(32*ECb))**(1/4)
phi0 = ((2*ECb)/ECj)**(1/4)
J = (EJ**2)/(2*ECj)
Final = []
   
k = np.arange(3,13)

def randn_complex(size):
    return (np.random.uniform(0,1, size) + 1j*np.random.uniform(0, 1, size)) / np.sqrt(2)

for n in k:

    qD1 = int((n*(n-1))/2-1)
    qD2 = 1
    qD3 = n - 1

    # for hamiltonian, not adding identities on either side and comparing the outcome to out function with mpo


    def F1(a1, a2, n):

        # now, calculating the E energy matrix
        EC0 = ECb
        EC1 = np.square(ECb)/(2*ECj)

        k = [EC1*np.ones(n-1),EC0*np.ones(n),EC1*np.ones(n-1)]
        offset = [-1,0,1]
        E = diags(k,offset).toarray()
        E[0, -1] = -EC1
        E[-1, 0] = -EC1
        # print(E)
        
        T = []
        for i in range(n):
            for j in range(n-i-1):
                if E[i,i+j+1] != 0:
                    t = [np.identity(2)]*(j+2)
                    t[0] = np.multiply(a1, (4*(n0**2)*E[i,i+j+1]))
                    t[-1] = np.multiply(a2, (4*(n0**2)*E[i,i+j+1]))
                    T.append(t)
                else:
                    continue
                        
        return T

    l_1 = [[a, aH], [-a, -aH], [aH, a], [-aH, -a], [-a, -aH], [a, aH], [-a, -aH], [aH, a]]


    def F2(a1, a2, n):

        T = [np.multiply(J,a1),np.multiply(J, a2)]
        return T

    a1 = np.tensordot(a, aH, axes = 1)
    a2 = np.tensordot(np.square(a), np.square(aH), axes = 1)
    a3 = np.tensordot(np.square(a), aH, axes = 1)
    a4 = np.tensordot(a, np.square(aH), axes = 1)


    l_2 = [[I, -a1], [aH, a], [a, aH], [-a1, I], np.multiply(1/4,[I, a2]),
            np.multiply(-1/2,[aH, a3]), np.multiply(-1/2,[a, a4]), np.multiply(1/4,[np.square(aH), np.square(a)]), [a1, a1],
            np.multiply(1/4,[np.square(a), np.square(aH)]), np.multiply(-1/2,[a4, a]),np.multiply(-1/2,[a3, aH]),
            np.multiply(-1/4,[a2, -I])]
    

    
    
    def F3(a1, a2, n):

        T = [np.identity(2)] * n
        T[0] = np.multiply(a1, J)
        T[-1] = np.multiply(a2, J)     
            
        return T

    l_3 = [[I, -a1], [-aH, -a], [-a, -aH], [-a1, I], np.multiply(1/4,[I, a2]),
            np.multiply(1/2,[aH, a3]), np.multiply(1/2,[a, a4]), np.multiply(1/4,[np.square(aH), np.square(a)]), [a1, a1],
            np.multiply(1/4,[np.square(a), np.square(aH)]), np.multiply(1/2,[a4, a]),np.multiply(1/2,[a3, aH]),
            np.multiply(1/4,[a2, I])]
   
   
    f1 = []
    for i in range(len(l_1)):
        site = 0
        x = F1(l_1[i][0], l_1[i][-1], n)
        y = F2(l_1[i][-1], l_1[i][0], n)
        z = F3(l_1[i][-1], l_1[i][0], n)

        f1.append(OpChain(x[0], [0]*(len(x[0])-1), istart = site))
        f1.append(OpChain(z, [0]*qD3, istart = site))
        for site in range(n-1):
            j = site + 1
            f1.append(OpChain(x[j], [0]*(len(x[j])-1), istart = site))
            f1.append(OpChain(y, [0], istart = site))


    f2 = []
    for i in range(len(l_2)):
        site = 0
        f2.append(OpChain(F3(np.multiply(l_2[i][-1],phi0**2), np.multiply(l_2[i][0],phi0**2), n), [0]*qD3, istart = site))
        for site in range(n-1):
            if i < 4:
                f2.append(OpChain(F2(np.multiply(l_2[i][0],phi0**2), np.multiply(l_2[i][-1],phi0**2), n), [0], istart = site))
                
            else:
                f2.append(OpChain(F2(np.multiply(l_2[i][0],phi0**4), np.multiply(l_2[i][-1],phi0**4), n), [0], istart = site))
                
                
    f3 = []
    for i in range(len(l_3)):

        if i < 4:
            f3.append(OpChain(F3(np.multiply(l_3[i][0],phi0**2), np.multiply(l_3[i][-1],phi0**2), n), [0]*qD3, istart = 0))
            
        else:
            f3.append(OpChain(F3(np.multiply(l_3[i][0],phi0**4), np.multiply(l_3[i][-1],phi0**4), n), [0]*qD3, istart = 0))

    F = []
    F.extend(f1)
    F.extend(f2)
    F.extend(f3)

    mpo = MPO.from_opchains([0, 0], n, F)

    D = [1] + [10]*(n-1) + [1]
    # create tensors with random complex values
    psi = MPS(mpo.qd, [np.zeros(Di, dtype=int) for Di in D], fill='postpone')
    np.random.seed(seed = 587439)
    random_3tens = randn_complex((2, 10, 10))
    for i in range(len(psi.A)):
        psi.A[i] = random_3tens
    psi.A[0] = psi.A[0][:, 0, :]
    psi.A[0] = np.reshape(psi.A[0] ,(2,1,10))
    psi.A[-1]= psi.A[-1][:, :, 0]
    psi.A[-1] = np.reshape(psi.A[-1], (2,10,1))

    numsweeps = 30
    final = minimization.calculate_ground_state_local_singlesite(mpo, psi, numsweeps, 25)
    Final.append(final[-1])


'''now, trying to plot the graph of frequency (en_min / h) against n (nr of capacitors)'''
# plotting the points
plt.plot(k, Final)
# plt.plot(range(numsweeps), final, 'b')
# plt.plot(range(D[i]), Final, 'r')

plt.xlabel('number of N capacitors')
plt.ylabel('minimum energy')

# function to show the plot
plt.show()