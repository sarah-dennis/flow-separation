# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import readwrite as rw
import stokes_pressure as pressure

def linf(err_ijs, N):
    return np.max(err_ijs)

def l1(err_ijs, N):
    return np.sum(err_ijs)

def l2(err_ijs, N):
    return np.sqrt(np.sum([e**2 for e in err_ijs]))


def linf_2D(err_ijs_u, err_ijs_v, N):
    return max(np.max(err_ijs_u), np.max(err_ijs_v))

def l1_2D(err_ijs_u, err_ijs_v, N):
    return (np.sum(err_ijs_u) + np.sum(err_ijs_v))/(2*N)

def l2_2D(err_ijs_u, err_ijs_v, N):
    return np.sqrt(np.sum([e_u**2 + e_v**2 for (e_u, e_v) in zip(err_ijs_u, err_ijs_v)])/(2*N))


def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates

#-------------------------------------------------------------------------------
def stokes_cnvg_self(Ex, args, U, Q, Re, N_min, Ns, N_max):
    ex_min = Ex(args, U, Q, Re, N_min)
    ex_max = Ex(args, U, Q, Re, N_max)
    size_min=ex_min.Ny*ex_min.Nx
    
    # Load max grid for 'true' values
    u_max, v_max, psi_max = rw.read_stokes(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))

    p_max = pressure.pressure(ex_max, u_max, v_max)
    dp_max = pressure.get_dp(ex_max, p_max) 
    p_max = p_max.reshape((ex_max.Ny,ex_max.Nx))

    u_max = u_max.reshape((ex_max.Ny, ex_max.Nx))
    v_max = v_max.reshape((ex_max.Ny, ex_max.Nx))

    mult_max = int(N_max/N_min)
    
    err_psi = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_psi_inf = np.zeros(len(Ns)+1)
    err_psi_l1 = np.zeros(len(Ns)+1)
    err_psi_l2 = np.zeros(len(Ns)+1)
    
    err_p = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_p_inf = np.zeros(len(Ns)+1)
    err_p_l1 = np.zeros(len(Ns)+1)
    err_p_l2 = np.zeros(len(Ns)+1)
    err_dp = np.zeros(len(Ns)+1)

    err_u = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_v = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_uv_inf = np.zeros(len(Ns)+1)
    err_uv_l1 = np.zeros(len(Ns)+1)
    err_uv_l2 = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = ex_min
            mult = 1
        else:
            ex_n = Ex(args, U, Q, Re, Ns[n-1])
            mult = int(Ns[n-1]/N_min)
            
            
        u_n, v_n, psi_n = rw.read_stokes(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
        psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))
        
        p_n = pressure.pressure(ex_n, u_n, v_n)
        dp_n = pressure.get_dp(ex_n, p_n)
        p_n=p_n.reshape((ex_n.Ny,ex_n.Nx))

        u_n = u_n.reshape((ex_n.Ny, ex_n.Nx))
        v_n = v_n.reshape((ex_n.Ny, ex_n.Nx))
        
        
        
        for k_min in range(size_min):
            # indices (i,j) on grid N_min
            i = k_min % ex_min.Nx
            j = (k_min // ex_min.Nx)

            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            
            if ex_max.space[j_max,i_max] != 1 or ex_n.space[j, i] !=1:
                continue
            
            err_psi[n,j,i] = abs(psi_max[j_max,i_max] - psi_n[j_n,i_n])
            err_p[n,j,i] = abs(p_max[j_max,i_max] - p_n[j_n,i_n])
            err_u[n,j,i]= abs(u_max[j_max,i_max]-u_n[j_n,i_n])
            err_v[n,j,i]= abs(v_max[j_max,i_max]-v_n[j_n,i_n])

    
        err_psi_inf[n] = linf(err_psi[n], size_min) 
        err_psi_l1[n] = l1(err_psi[n], size_min)
        err_psi_l2[n] = l2(err_psi[n], size_min)
        
        err_p_inf[n] = linf(err_p[n], size_min) 
        err_p_l1[n] = l1(err_p[n], size_min) 
        err_p_l2[n] = l2(err_p[n], size_min)
        err_dp[n] = abs(dp_max-dp_n)
        
        err_uv_inf[n] = linf_2D(err_u[n], err_v[n], size_min) 
        err_uv_l1[n] = l1_2D(err_u[n], err_v[n], size_min)    
        err_uv_l2[n] = l2_2D(err_u[n], err_v[n], size_min)    
            
            
    l1_rate = convg_rate(err_psi_l1)
    l2_rate = convg_rate(err_psi_l2)
    inf_rate = convg_rate(err_psi_inf)
    cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)
    
   
    p_l1_rate = convg_rate(err_p_l1)
    p_l2_rate = convg_rate(err_p_l2)
    p_inf_rate = convg_rate(err_p_inf)
    dp_rate = convg_rate(err_dp)
    p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, dp_rate], axis=0)
    

    uv_l1_rate = convg_rate(err_uv_l1)
    uv_l2_rate = convg_rate(err_uv_l2)
    uv_inf_rate = convg_rate(err_uv_inf)
    uv_cnvg_rates = np.stack([uv_l1_rate, uv_l2_rate, uv_inf_rate], axis=0)
    
        
    print("stream cnvg rates")
    print("l1: " + np.array2string(cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(cnvg_rates[1], precision=2))
    print("linf" + np.array2string(cnvg_rates[2], precision=2))

    print("pressure cnvg rates")
    print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
    print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
    print("dp: " + np.array2string(p_cnvg_rates[3], precision=2))
    

    print("velocity cnvg rates")
    print("l1: " + np.array2string(uv_cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(uv_cnvg_rates[1], precision=2))
    print("linf" + np.array2string(uv_cnvg_rates[2], precision=2))
    
        
    stream_errs=[err_psi_l1, err_psi_l2, err_psi_inf]
    p_errs = [err_p_l1, err_p_l2, err_p_inf]
    uv_errs = [err_uv_l1, err_uv_l2, err_uv_inf]
    return stream_errs, p_errs, uv_errs


# def stokes_cnvg_self(Ex, args, U, Q, Re, N_min, Ns, N_max):
#     ex_min = Ex(args, U, Q, Re, N_min)
#     ex_max = Ex(args, U, Q, Re, N_max)
#     # size_min=ex_min.Ny*ex_min.Nx
    
#     # Load max grid for 'true' values
#     u_max, v_max, psi_max = rw.read_stokes(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
#     psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))

#     p_max = pressure.pressure(ex_max, u_max, v_max)
#     dp_max = pressure.get_dp(ex_max, p_max) 
#     p_max = p_max.reshape((ex_max.Ny,ex_max.Nx))

#     u_max = u_max.reshape((ex_max.Ny, ex_max.Nx))
#     v_max = v_max.reshape((ex_max.Ny, ex_max.Nx))


#     # mult_max = int(N_max/N_min)
    
#     err_psi_inf = np.zeros(len(Ns)+1)
#     err_psi_l1 = np.zeros(len(Ns)+1)
#     err_psi_l2 = np.zeros(len(Ns)+1)
#     err_p_inf = np.zeros(len(Ns)+1)
#     err_p_l1 = np.zeros(len(Ns)+1)
#     err_p_l2 = np.zeros(len(Ns)+1)
#     err_dp = np.zeros(len(Ns)+1)
#     err_uv_inf = np.zeros(len(Ns)+1)
#     err_uv_l1 = np.zeros(len(Ns)+1)
#     err_uv_l2 = np.zeros(len(Ns)+1)
    
#     for n in range(len(Ns)+1):
        
#         # Load grid at N < N_max
#         if n == 0:
#             ex_n = ex_min
#             scale = int(N_max/N_min)
#         else:
#             ex_n = Ex(args, U, Q, Re, Ns[n-1])
#             scale = int(N_max/Ns[n-1])

#         size_n = ex_n.Nx * ex_n.Ny
#         u_n, v_n, psi_n = rw.read_stokes(ex_n.filestr+".csv", size_n)
#         psi_n = psi_n.reshape((ex_n.Ny, ex_n.Nx))

#         p_n = pressure.pressure(ex_n, u_n, v_n)
#         dp_n = pressure.get_dp(ex_n, p_n)
#         p_n=p_n.reshape((ex_n.Ny,ex_n.Nx))

#         u_n = u_n.reshape((ex_n.Ny, ex_n.Nx))
#         v_n = v_n.reshape((ex_n.Ny, ex_n.Nx))
        
#         err_psi_n = np.zeros((ex_n.Ny,ex_n.Nx))
#         err_u_n = np.zeros((ex_n.Ny,ex_n.Nx))
#         err_v_n = np.zeros((ex_n.Ny,ex_n.Nx))
#         err_p_n = np.zeros((ex_n.Ny,ex_n.Nx))

       
#         for k_n in range(size_n):
            
#             # indices on grid n
#             i = k_n % ex_n.Nx
#             j = k_n // ex_n.Nx

#             # -> indices on grid N_max || N_min
#             i_max = scale * i
#             j_max = scale * j
            
#             if ex_max.space[j_max,i_max] != 1 or ex_n.space[j, i] !=1:
#                 continue
            
#             err_psi_n[j,i] = abs(psi_max[j_max,i_max] - psi_n[j,i])
#             err_p_n[j,i] = abs(p_max[j_max,i_max] - p_n[j,i])
#             err_u_n[j,i] = abs(u_max[j_max,i_max] - u_n[j,i])
#             err_v_n[j,i] = abs(v_max[j_max,i_max] - v_n[j,i])

    
#         err_psi_inf[n] = linf(err_psi_n, size_n) 
#         err_psi_l1[n] = l1(err_psi_n, size_n)
#         err_psi_l2[n] = l2(err_psi_n, size_n)
        
#         err_p_inf[n] = linf(err_p_n, size_n) 
#         err_p_l1[n] = l1(err_p_n, size_n) 
#         err_p_l2[n] = l2(err_p_n, size_n)
#         err_dp[n] = abs(dp_max-dp_n)
        
#         err_uv_inf[n] = linf_2D(err_u_n, err_v_n, size_n) 
#         err_uv_l1[n] = l1_2D(err_u_n, err_v_n, size_n)    
#         err_uv_l2[n] = l2_2D(err_u_n, err_v_n, size_n)    
            
            
#     l1_rate = convg_rate(err_psi_l1)
#     l2_rate = convg_rate(err_psi_l2)
#     inf_rate = convg_rate(err_psi_inf)
#     cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)
    
   
#     p_l1_rate = convg_rate(err_p_l1)
#     p_l2_rate = convg_rate(err_p_l2)
#     p_inf_rate = convg_rate(err_p_inf)
#     dp_rate = convg_rate(err_dp)
#     p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, dp_rate], axis=0)
    

#     uv_l1_rate = convg_rate(err_uv_l1)
#     uv_l2_rate = convg_rate(err_uv_l2)
#     uv_inf_rate = convg_rate(err_uv_inf)
#     uv_cnvg_rates = np.stack([uv_l1_rate, uv_l2_rate, uv_inf_rate], axis=0)
    
        
#     print("stream cnvg rates")
#     print("l1: " + np.array2string(cnvg_rates[0], precision=2))
#     print("l2: " + np.array2string(cnvg_rates[1], precision=2))
#     print("linf" + np.array2string(cnvg_rates[2], precision=2))

#     print("pressure cnvg rates")
#     print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
#     print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
#     print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
#     print("dp: " + np.array2string(p_cnvg_rates[3], precision=2))
    

#     print("velocity cnvg rates")
#     print("l1: " + np.array2string(uv_cnvg_rates[0], precision=2))
#     print("l2: " + np.array2string(uv_cnvg_rates[1], precision=2))
#     print("linf" + np.array2string(uv_cnvg_rates[2], precision=2))
    
        
#     stream_errs=[err_psi_l1, err_psi_l2, err_psi_inf]
#     p_errs = [err_p_l1, err_p_l2, err_p_inf]
#     uv_errs = [err_uv_l1, err_uv_l2, err_uv_inf]
#     return stream_errs, p_errs, uv_errs



#------------------------------------------------------------------------------
    

#------------------------------------------------------------------------------


