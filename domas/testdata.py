import numpy as np
import matplotlib.pyplot as plt

n = 50

# 0         1           2       3       4       5           6       
#alpha[deg] FFl[lbhr] FFr[lbhr] hp[ft] mfu[lb] Tmta[c] Vc[kts] 

i = [700,860,
     1780,1890,
     3030,3200,
     5460,5500,
     7140,7340,
     7860,7980]

a = np.genfromtxt('datamatlab.txt')

plt.subplot(311)
for j in i:
    plt.axvline(j)
plt.plot(a[3])

plt.subplot(312)
for j in i:
    plt.axvline(j)
plt.plot(a[6])

plt.subplot(313)
for j in i:
    plt.axvline(j)
for j in [ 1.3,  2. ,  3. ,  4.5,  8.3, 10.1]:
    plt.axhline(j)
plt.plot(a[0])

print('Close plot to write.')
plt.show()

#for i in [ 1.3,  2. ,  3. ,  4.5,  8.3, 10.1]:
#    plt.axhline(i)

t = a.T


y = np.array(t[i[0]:i[1]])


# From velocity plot, read flat
for j in range(2,len(i),2):
    y = np.concatenate((y, t[i[j]:i[j+1]]))
print(len(y),'data points')
a = y.T

np.savetxt('dataflat.txt',a)
print('Written to file.')

