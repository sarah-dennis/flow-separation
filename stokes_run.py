# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples


zoom_on = True   #set zoom location in graphics.py

#------------------------------------------------------------------------------

U=0
Q=1


# U=1
# Q=0

Re=0


#------------------------------------------------------------------------------
h_in = 1
h_out =2 #2.75, 2, 1.25
l_in = 8
l_out = 8
args = [h_in, h_out, l_in, l_out]
Example = examples.BFS


#------------------------------------------------------------------------------
# h_in = 1
# h_out = 2
# l_in = 8
# l_out=8

# xr = 0.175 #  0.35, 0.2625, 0.175
# yr = 0.2  #  0.4,  0.3,    0.2

# args = [h_in, h_out, l_in, l_out, xr, yr]
# Example = examples.BFS_wedge

#------------------------------------------------------------------------------

# hin=1
# hout=2
# L=16
# delta=1/4

# args = [hin, hout, L, delta]
# Example = examples.BFS_pwl


#------------------------------------------------------------------------------
# H = 1/2
# L = 1

# args = [H, L]# tri slope = 2H/L
# Example = examples.TriCavity

#------------------------------------------------------------------------------
solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=50000)                

N=160
# solver.new_run(N) 
# solver.load_run(N)

# solver.load_scale(N,2*N) 
# solver.load_run(2*N)
# solver.load_copy(N, new_Example, new_args)

# solver.load_run_many(N, 2, 5)

# solver.new_run_many(N, 2, 4)  
# solver.load_run_new_many(N, 2,2)

solver.load_plot(N, zoom=zoom_on)

# ------------------------------------------------------------------------------
solver.compare(args, U, Q, Re, 20, [40,80,160],320, p_err= False, uv_err=True)







