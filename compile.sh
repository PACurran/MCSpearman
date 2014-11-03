#!/bin/bash
# 
# Requires mcspearman_v#.c, mcspearman.h 
#
# Requires the GNU Scientific Library (GSL) development & library packages


gcc -Wall -o  mcspearman mcspearman_v0.3.c -L/usr/lib -lgsl -lgslcblas -lm

# where -L/usr/lib is your GSL directory 
