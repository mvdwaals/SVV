import numpy as np
import matplotlib.pyplot as plt

n = 50

a = np.genfromtxt('datamatlab.txt')
t = a.T

t1 = t[0:860]
t2 = t[1780:2040]
t3 = t[2800:2940]
t4 = t[3800:4040]
t5 = t[5280:5380]
t6 = t[7750:8300]

t1 = t1[0::int(170/n)+1]
t2 = t2[0::int(152/n)+1]
t3 = t3[0::int(128/n)+1]
t4 = t4[0::int(148/n)+1]
t5 = t5[0::int(120/n)+1]
t6 = t6[0::int(1110/n)+1]

y = np.concatenate((t1, t2, t3, t4, t5, t6))
print(len(y))
a = y.T

np.savetxt('dataflat.txt',a)

# 0         1           2       3       4       5           6       
#alpha[deg] FFl[lbhr] FFr[lbhr] hp[ft] mfu[lb] Tmta[c] Vc[kts] 

