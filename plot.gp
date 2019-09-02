#!/usr/bin/env gnuplot
set terminal pdfcairo fontscale 0.7


set pm3d
set pm3d corners2color c1
unset surface
set view map

set rmargin at screen 0.8
set lmargin at screen 0.15

set grid
set xlabel "x (mm)"
set ylabel "y (mm)"
set cblabel "<B,n_{sample}> (T)" 
set output 'data_tilt_x=1_sample_plain.pdf'
B0 = 0.1065
splot 'data_tilt_x=1_sample_plain.dat' using (1000 * $1):(1000 * $2):($10-B0)/B0*100 title ''

set output 'data_tilt_x=0_sample_plain.pdf'
splot 'data_tilt_x=0_sample_plain.dat' using (1000 * $1):(1000 * $2):($10-B0)/B0*100 title ''
