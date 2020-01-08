#!/usr/bin/env gnuplot
set terminal pdfcairo fontscale 0.7


set pm3d
set pm3d corners2color c1
unset surface
set view map

set rmargin at screen 0.75
set lmargin at screen 0.15

set grid
set xlabel "x (mm)"
set ylabel "y (mm)"
set cblabel "((b-b_0)/b_0)*100; b = B•n_{sample}" 

set output 'data_tilt_x=1_sample_plain.pdf'

B0 = 0.1065 # field in (0,0,0)

splot 'data_tilt_x=1_sample_plain.dat' using (1000 * $1):(1000 * $2):($10-B0)/B0*100 title ''

set output 'data_tilt_x=0_sample_plain.pdf'
splot 'data_tilt_x=0_sample_plain.dat' using (1000 * $1):(1000 * $2):($10-B0)/B0*100 title ''


set output 'data_tilt_x=0_shift_z=0.002_sample_plain.pdf'
splot 'data_tilt_x=0_shift_z=0.002_sample_plain.dat' using (1000 * $1):(1000*$2):($10-B0)/B0*100 title ''