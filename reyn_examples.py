# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWC_Height,  PWL_Height

# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        |______|
#


class BFS(PWC_Height):
    def __init__(self, args, N):
        h_in, h_out, l_in, l_out = args
        x0 = 0
        xf = l_in+l_out
        N_regions = 2
        x_peaks = np.asarray([0, l_in, xf], float)
        h_peaks = np.asarray([[h_in, h_in], [h_in, h_out], [h_out, h_out]], float)
        namestr = ''#f'BFS_H{int(H)}L{int(xf)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)


# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        \______|
#

        
class BFS_pwl(PWL_Height):
    def __init__(self, args, N):
        H, h, L, delta = args

        x0 = 0
        xf = L
        N_regions = 3
        x_peaks = np.asarray([x0, (L-delta)/2, (L+delta)/2, xf], float)
        h_peaks = np.asarray(
            [[0, H], [H, H], [h,h], [h, 0]], float)
        namestr = ''#f'dBFS_H{int(H)}L{int(xf)}_d{int(delta)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)


class BFS_wedge(PWL_Height):
    def __init__(self, args, N):
        H, h, la, lb, xr, yr = args
        x0 = 0
        xf = la+lb
        N_regions = 3
        x_peaks = np.asarray([x0, la-xr, la, xf], float)
        h_peaks = np.asarray( [[0, H], [H, H], [H-yr, h], [h, 0]], float)
        namestr = ''#f'cBFS_H{int(H)}L{int(xf)}_xr{int(xr)}yr{int(yr)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)

# -----------------------------------------------------------------------------------------------------------------------------------
#   ______
#   \    /
#    \  /
#     \/


class TriCavity(PWL_Height):
    def __init__(self, args, N):
        H, L = args
        l_a = L/2
        l_b = L/2
        x0 = 0
        xf = x0 + l_a + l_b
        N_regions = 2

        H = args[0]

        x_peaks = np.asarray([x0, x0 + l_a, xf], float)

        dx = 0.01
        # h = max(args[2],dx) #H/2  # dx
        
        h_peaks = np.asarray(([[0, dx], [H, H], [dx, 0]]), float)
        namestr ='' #'TriCavity_H{int(H)}L{int(xf)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)
 