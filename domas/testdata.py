import numpy as np
import matplotlib.pyplot as plt

n = 50

a = np.genfromtxt('datamatlab.txt')

b = np.gradient(a[0])
c = np.gradient(a[3])
d = np.gradient(a[6])

#plt.plot(a[0])
plt.plot(np.sqrt(c**2 + d**2) < 0.1)

#for i in [ 1.3,  2. ,  3. ,  4.5,  8.3, 10.1]:
#    plt.axhline(i)

t = a.T

i = [700,860,
     1780,2040,
     3030,3200,
     5460,5500,
     7140,7340,
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

