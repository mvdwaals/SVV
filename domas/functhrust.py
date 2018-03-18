from dataconst import *
import os
import numpy as np

def fTtotal(T, n_tests, hp, M, FFl, FFr):
    TISA = T0 + labda * np.array(hp)
    DeltaT = T - TISA
    f = open('matlab.dat','w')
    for i in range(n_tests):
        str1 = str(hp[i])+'\t'
        str2 = str(M[i])+'\t'
        str3 = str(DeltaT[i])+'\t'
        str4 = str(FFl[i])+'\t'
        str5 = str(FFr[i])+'\r\n'
        strtot = str1+str2+str3+str4+str5
        f.write(strtot)
    f.close()

    os.system('java -jar thrust.jar')

    gen = np.genfromtxt('thrust.dat').T
    T1 = gen[0]
    T2 = gen[1]
    return T1 + T2
