# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""


import reyn_boundary as bc
import reyn_examples as examples
import reyn_solvers as solvers

#-------------------plotting---------------------------------------------------

plots_on = True
uv_on = False          # plot u(x,y) & v(x,y) & |(u,v)|
inc_on = False         # plot ux + vy =? 0
zoom_on =False         # plot a zoomed-in window, set location in reyn_solution.py

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

# Example = examples.BFS
# H=1 
# h=2
# l=8
# l_out=8
# args =  [h, H, l, l_out]


Example = examples.BFS_pwl
H = 2
h=1
delta = 1
L=16
args = [H,h,L,delta]


# Example = examples.BFS_wedge
# h = 1
# H = 2
# l = 3
# L = 3
# xr = 0.5
# yr = 0.5
# args = [H, h, l, L, xr, yr]


# Example = examples.TriCavity
# H=4 # apex height
# L=2
# args = [H, L]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

#fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 8
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = solvers.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 100
# solution = solver.fd_solve(N)

solution = solver.pwl_solve(N)


if plots_on:
    solution.p_plot(zoom=zoom_on)
    solution.v_plot(zoom=zoom_on, uv=uv_on, inc=inc_on)
#------------------------------------------------------------------------------
