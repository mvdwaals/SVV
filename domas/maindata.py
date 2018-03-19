import numpy as np
import matplotlib.pyplot as plt

n = 50

a = np.genfromtxt('datamatlab.txt')
plt.plot(a[3])
plt.plot(a[6]+5700)

t = a.T

i = [0,860,
     1780,2040,
     3250,3330,
     3770,3860,
     5460,5550,
     7860,7980]

for j in i:
    plt.axvline(j)
plt.show()

# From velocity plot, read flat
t1 = t[i[0]:i[1]]
t2 = t[i[2]:i[3]]
t3 = t[i[4]:i[5]]
t4 = t[i[6]:i[7]]
t5 = t[i[8]:i[9]]
t6 = t[i[10]:i[11]]


y = np.concatenate((t1, t2, t3, t4, t5, t6))
#y = t # Uncomment this and run to use all matlab data (not recommended)
print(len(y))
a = y.T

np.savetxt('dataflat.txt',a)
print('Written to file.')

# 0         1           2       3       4       5           6       
#alpha[deg] FFl[lbhr] FFr[lbhr] hp[ft] mfu[lb] Tmta[c] Vc[kts] 

