#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import loopfield


# compensation coils:
# r_in = 5mm, r_out = 30.25mm
# d_coil = 5.3mm
# w = 25mm + d_coil/2

# N_loop_1 = 9530
# N_loop_2 = 10623
# d_nomm = 100μm

# coil axis in z-direction
# sample midpoint in (0,0,0)
# untilted sample normal to (0,0,1)


field = loopfield.Field()


normal = [0,0,1.]

def add_loop(position, normal, radius, current):
    print("adding loop at ", position, ", with radius ", radius)
    field.addLoop(loopfield.Loop(position, normal, radius, current))


    
def add_coil(pos, N, normal=[0,0,1], width=5.3e-3, r_o=30.25e-3, r_i=5e-3, current=1):
    # pos is coil midpoint
    A = width * (r_o - r_i)
    A_loop = A/N
    d_loop = np.sqrt(A_loop)
    N_d = int(round(width/d_loop))
    N_r = int(round((r_o-r_i)/d_loop))
    print("building coil with %d x %d = %d loops, error = %.2g" % (N_d, N_r, N_d * N_r, np.abs((N_d*N_r - N) / N)))
    radiuses = np.linspace(r_i, r_o, N_r)
    mid_points = pos + np.outer(np.linspace(-width/2, width/2, N_d), normal)
    print(mid_points)
    for radius in (radiuses):
        for mid_point in (mid_points):
            add_loop(mid_point, normal, radius, current)
    
coil_pos = np.array([0, 0, 25e-3])
add_coil(coil_pos, 9530)
add_coil(-coil_pos, 10623)


# # calculate field on axis
# N = 100
# pos_vals = np.zeros((N,3))
# z_vals = np.linspace(-1.5*w, 1.5*w, N)
# pos_vals[:,2] = z_vals
# field_vals = field.evaluate(pos_vals)
# print(pos_vals)
# plt.plot(z_vals, field_vals[:,2])
# plt.grid()
# plt.show(block=False)

import code
code.interact()
