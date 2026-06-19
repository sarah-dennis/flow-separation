# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_examples as examples
import graphics  
import readwrite as rw



def get_bfs_attachments(ex):
       
    u, v, psi = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    
    y_xr = ex.yf   # xr on y=h
    x_yr = ex.xf/2 # yr on x=L/2
    
    j_xr =int(y_xr/ex.dy)-1
    i_yr =int(x_yr/ex.dx)+1 
    

    xrs = []
    yrs = []

    # find xr
    for i in range(ex.Nx-1):
        k_a = int(j_xr*ex.Nx + i)
        k_b = int(j_xr*ex.Nx + i+1 )   
        if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):

            xr = ex.xs[i]-x_yr
            if xr > 0:
                xrs.append(xr)
    #find yr        
    for j in range(ex.Ny-1):
        k_a = int(j*ex.Nx + i_yr)
        k_b = int((j+1)*ex.Nx + i_yr)
        
        if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
            
            yr = y_xr-ex.ys[j] 
            if yr > 0:
                yrs.append(yr)
                
    if xrs==[]:
        xrs.append(0)
    if yrs == []:
        yrs.append(0)
        
    return xrs, yrs


def get_bfs_separation(ex):
       
    u, v, psi = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)

    xs = []
    ys = []

    # find xr
    for i in range(ex.Nx//2+1, ex.Nx-1):
        for j in range(ex.Ny-5):

            k_a = int(j*ex.Nx + i)
            k_b = int((j+1)*ex.Nx + i)   
            
            if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
                xs.append(ex.xs[i])
                ys.append(ex.ys[j])
                continue #primary sep only
           
   
    return xs, ys


U=0
Q=1
Re=0
#------------------------------------------------------------------------------

Example = examples.BFS
h_outs = [1.25, 2, 2.75] #N=160

h_in = 1
l_in = 8
l_out = 8

args_all = [[h_in, h_out, l_in, l_out] for h_out in h_outs]

num_tests = len(h_outs)

#------------------------------------------------------------------------------
# xrs = np.zeros(num_tests)
# yrs = np.zeros(num_tests)
# N = 160

# for k in range(num_tests):
#     args = args_all[k]         
#     ex = Example(args, U, Q, Re, N)
#     xr, yr = get_bfs_attachments(ex)
    
#     print(h_outs[k], xr, yr)
#     #primary
#     xrs[k] = xr[-1]
#     yrs[k] = yr[0]

# # print(xrs,yrs)

# graphics.plot_2D_multi([xrs, yrs], h_outs, 'BFS Points of Flow Separation', ['$x_r$', '$y_r$'], ['$\mathcal{H}$', 'length'], loc='left')
#------------------------------------------------------------------------------
N = 160

for k in range(num_tests):
    args = args_all[k]         
    ex = Example(args, U, Q, Re, N)
    xrs, yrs = get_bfs_attachments(ex)
    xr = xrs[-1]
    yr = yrs[0]
    xs, ys = get_bfs_separation(ex)
    lin_ys = [yr/xr * (x-(ex.xf/2))+ (ex.yf-yr)-ex.dy for x in xs]
    
    n = len(xs)
    errs = [abs(y - lin_y) for (y, lin_y) in zip(ys, lin_ys)]
    mae = np.sum(errs)/n
    sqrs = [err**2 for err in errs]
    rmse = np.sqrt(np.sum(sqrs)/n)
    
    
    print('H: %.4f, xr: %.4f, MAE: %.4f, RMSE: %.4f'%(h_outs[k], xr, mae, rmse))

    graphics.plot_2D_multi([ys, lin_ys], xs, 'flow partition', ['true', 'est'], ['$x$', '$y$'], loc='right')

