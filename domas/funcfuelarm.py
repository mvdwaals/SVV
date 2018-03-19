import numpy as np
from dataconst import *

mass = [100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100., 2200., 2300., 2400., 2500., 2600., 2700.]

momentarm = [298.16, 591.18, 879.08, 1165.42, 1448.40, 1732.53, 2014.80, 2298.84, 2581.92, 2866.30, 3150.18, 3434.52, 3718.52, 4003.23, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.90, 5994.04, 6278.47, 6562.82, 6846.96, 7131.00, 7415.33, 7699.60]

def findarm(m):

    for i in range(len(mass)-1):
       if m >= mass[i] and m < mass[i+1]:
           diff = momentarm[i+1] - momentarm[i]
           slope = diff / (mass[i+1]-mass[i])
           x = (m - mass[i]) / (mass[i+1]-mass[i])
           momentarmnew = (momentarm[i] + x * diff)*inlbf_kgm/g/100
    return momentarmnew

