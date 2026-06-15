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
zoom_on = not False        # plot a zoomed-in window, set location in graphics.py

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

# Example = examples.BFS
# h_in=2
# h_out=1
# l_in=8
# l_out=8
# args =  [h_in, h_out, l_in, l_out]


# Example = examples.BFS_pwl
# h_in = 2
# h_out=1
# delta = 1/8
# L=16
# args = [h_in,h_out,L,delta]


Example = examples.BFS_wedge
h_in = 1
h_out = 2
l_in = 8
l_out = 8
xw = 0.35
yw = 0.4
args = [h_in, h_out, l_in, l_out, xw, yw]


# Example = examples.TriCavity
# H=4 
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


N = 80
# solution = solver.fd_lu_solve(N)

solution = solver.pwl_solve(N)

print('solve time: %.2fs'%solution.time)
if plots_on:
    solution.p_plot(zoom=zoom_on)
    solution.v_plot(zoom=zoom_on, uv=uv_on, inc=inc_on)
#------------------------------------------------------------------------------
