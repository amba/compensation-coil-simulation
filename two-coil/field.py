#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os.path
import argparse
import sys
from scipy.spatial.transform import Rotation

import loopfield

# two-coil setup
#  width=2e-3


field = loopfield.Field()


def save_3d_file(output_file, data, header):
    fh = open(output_file, 'w')
    fh.write(header + "\n")
    shape = data.shape
    for i in range(shape[0]):
        block = data[i]
        np.savetxt(fh, block, fmt="%.17g", delimiter="\t")
        fh.write("\n")
    fh.close()

    
def add_coil(pos, N, normal=[0,0,1], width=5.3e-3, r_o=30.25e-3 - 1e-3, r_i=5e-3, current=1):
    # pos is coil midpoint
    A = width * (r_o - r_i)
    A_loop = A/N
    d_loop = np.sqrt(A_loop)
    N_d = int(round(width/d_loop))
    N_r = int(round((r_o-r_i)/d_loop))
    N_prod = N_d * N_r
    print("building coil with %d x %d = %d loops, error = %d ( %.2g percent)" % (N_d, N_r, N_prod, N_prod - N, 100 * np.abs((N_prod - N) / N)))
    radiuses = np.linspace(r_i, r_o, N_r)
    mid_points = pos + np.outer(np.linspace(-width/2, width/2, N_d), normal)
    for radius in (radiuses):
        for mid_point in (mid_points):
            field.addLoop(loopfield.Loop(mid_point, normal, radius, current))


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")
parser.add_argument('-o', '--output', help="basename for output data files")
parser.add_argument('--tilt_x', help="tilt (in degrees) around x axis", type=float, default=0)
parser.add_argument('--shift_z', help="shift sample in z direction (m)", type=float, default=0) 


args = parser.parse_args()

if args.output:
    args.output += "_tilt_x=%.2g" % args.tilt_x
    args.output += "_shift_z=%.2g" % args.shift_z
    
w = 25e-3 + 2e-3
coil_pos = np.array([0, 0, w])
add_coil(coil_pos, 9530)
add_coil(-coil_pos, 10623)

# calculate field on axis
N = 100
pos_vals = np.zeros((N,3))
z_vals = np.linspace(-2e-3, 2e-3, N)
pos_vals[:,2] = z_vals
field_vals = field.evaluate(pos_vals)
data_axis = np.stack((z_vals, field_vals[:,2]), axis=1)



# calculate field on sample plain
N = 50
d_sample = 6e-3
sample_normal = np.array([0,0,1])
x_vals = np.linspace(-d_sample/2, d_sample/2, N)
y_vals = np.linspace(-d_sample/2, d_sample/2, N)
x_mesh, y_mesh = np.meshgrid(x_vals, y_vals, indexing="ij")
sample_positions = np.zeros((N, N, 3))
sample_positions[...,0] = x_mesh
sample_positions[...,1] = y_mesh

# tilt sample

r = Rotation.from_euler('x', args.tilt_x, degrees=True)
rotation_matrix = r.as_dcm()

tilted_positions = np.dot(sample_positions, rotation_matrix.T)

# shift sample

tilted_positions = tilted_positions - np.array([0, 0, args.shift_z])

# calculate field

sample_fields = field.evaluate(np.reshape(tilted_positions, (N * N, 3)))
sample_fields = np.reshape(sample_fields, (N, N, 3))

sample_normal = np.dot(rotation_matrix, sample_normal)
print("tilted sample_normal = ", sample_normal)

sample_normal_field = np.dot(sample_fields, sample_normal)
sample_normal_field = np.reshape(sample_normal_field, (N, N, 1))
print(sample_normal_field.shape)

sample_data = np.concatenate((sample_positions, tilted_positions, sample_fields, sample_normal_field), axis=2)

def ensure_unique(file):
    if not args.force and os.path.isfile(file):
        sys.exit("file %s already exists. Use -f option to overwrite" % file)
        
if args.output:
    axis_output_file = args.output + "_axis.dat"
    ensure_unique(axis_output_file)
    header = "# z\tB_z"
    np.savetxt(axis_output_file, data_axis, header=header, comments='', fmt="%.17g")

    sample_output_file = args.output + "_sample_plain.dat"
    ensure_unique(sample_output_file)
    
    header = "# x\ty\tz\tx_t\ty_t\tz_t\tB_x\tB_y\tB_z\tB_normal"
    save_3d_file(sample_output_file, sample_data, header)




