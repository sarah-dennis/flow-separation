# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:10:36 2026

@author: sarah
"""
from time import time
import domain as dm
import reyn_boundary as bc
import reyn_velocity as rv
import reyn_pressure as rp
from reyn_control import Reyn_Solution
from reyn_heights import PWL_Height, PWC_Height, make_PWL


class Reynolds_Solver: 
    def __init__(self, Example, BC, args=None):
        self.Example = Example #initialize height = Example(args) in the solver
        self.args = args
        
        self.BC = BC        
        
#----------------------------------------------------------------------------------
    def fd_lu_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        t0 = time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        pressure = rp.FD_LU_ReynPressure(height, self.BC)
        tf = time()
            
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        t = tf-t0
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
#----------------------------------------------------------------------------------
#---------------------------------Piecewise----------------------------------------
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------   
    def pwl_solve(self, N):
        solver_title = "Reynolds" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        if not (isinstance(self.BC, bc.Mixed)): #TODO
            raise TypeError('Only prescribed flux for pwl schur solver')
        
        if not isinstance(height, PWL_Height) and not isinstance(height, PWC_Height):  # PWC is a PWL
            height = make_PWL(height)

        t0 = time()
        pressure = rp.PwlSchur_ReynPressure(height, self.BC)
        tf = time()
        t = tf-t0
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution  