# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_examples as examples
import graphics  
import readwrite as rw




def get_tri_centers(ex, H, L):
       
    u, v, psi = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    y_cen=[]
    psi_cen = []
    ic = int((L/2)/ex.dx)
    for j in range(1, ex.Ny-1):
        k_a = (j-1)*ex.Nx + ic
        k_b = j*ex.Nx + ic
        k_c = (j+1)*ex.Nx + ic
        # psi_cs[j] = psi[k]
        psi_a = psi[k_a]
        psi_b = psi[k_b]
        psi_c = psi[k_c]
        if psi_b > 0:
            if psi_b > psi_a and psi_b > psi_c:
                y_cen.append(j*ex.dx)
                psi_cen.append(psi_b)
        elif psi_b < 0:
            if psi_b < psi_a and psi_b < psi_c:
                y_cen.append(j*ex.dx)
                psi_cen.append(psi_b)
        
    return psi_cen, y_cen

def get_tri_attachments(ex, H, L):
       
    u, v, psi = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    xrs=[]
    
    for i in range(ex.Nx):
        if i < ex.Nx/2 - 1:
            j_ha = H/(L/2) * i
            j_hb = H/(L/2) * (i+1)
            if j_ha < ex.Ny and j_hb < ex.Ny:
                k_a = int(j_ha*ex.Nx + i)
                k_b = int(j_hb*ex.Nx + i+1)   
                
              
                if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
                    
                    xr = ex.xs[i]
                    xrs.append(xr)


    if xrs==[]:
        xrs.append(0)
        
    return xrs


def convg_tri_attachments(H, L, Ns):
    Example = examples.TriCavity
    args = [H, L]
    
    xr1s = []
    xr2s = []
    for N in Ns:
        ex = Example(args, U, Q, Re, N)
        xrs = get_tri_attachments(ex, H, L)
        xr1s.append(xrs[1])
        xr2s.append(xrs[2])

    err_xr1s = []
    xr1_best = xr1s[-1]
    for xr1 in xr1s[:-1]:
        err_xr1s.append(abs(xr1-xr1_best))
    
    err_xr2s = []
    xr2_best = xr1s[-1]
    for xr2 in xr2s[:-1]:
        err_xr2s.append(abs(xr2-xr2_best))
        
    # graphics.plot_2D_multi([err_xr1s, err_xr2s], Ns[:-1], '$x_r$ convergence',  ['$x_{r_1}$', '$x_{r_2}$'],['$N$', '$|x_{r^*}-x_r|$'])
    graphics.plot_2D(err_xr1s, Ns[:-1], '$x_r$ convergence', ['$N$', '$|x_{r^*}-x_r|$'])

U=1
Q=0
Re=0
#------------------------------------------------------------------------------


# hs = [0.25, 0.5, 0.75, 1, 1.5, 2, 3]
hs=[2]
L = 1
N=320
#------------------------------------------------------------------------------

num_h = len(hs)

for i in range(num_h):
    H= hs[i]

    Example = examples.TriCavity
    args = [H, L]
    
    
    ex = Example(args, U, Q, Re, N)
    
    xr = get_tri_attachments(ex, H, L)
    yr = [H/(L/2) * x for x in xr]
    print(' H ', H, '\n xr ', xr, '\n yr ', yr)

    psi_cen, y_cen = get_tri_centers(ex, H, L)
    print('psi ', psi_cen, '\n y ', y_cen)

#------------------------------------------------------------------------------

# H = 2
# L=1
# Ns = [20, 40, 80, 160]
# convg_tri_attachments(H, L, Ns)
    