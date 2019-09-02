#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

N = 9530
d = 5.3e-3
r_i = 5e-3
r_o = 30.25e-3
A = d * (r_o - r_i)
A_loop = A/N
d_loop = np.sqrt(A_loop)

print(A_loop, d_loop)

N_d = int(round(d/d_loop))
N_r = int(round((r_o-r_i)/d_loop))
print(N_d)
print(N_r)
print(N_d * N_r)

