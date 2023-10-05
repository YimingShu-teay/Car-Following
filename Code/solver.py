#MPC Model
import numpy as np
from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import math
from math import factorial
from celluloid import Camera
from pyomo.util.infeasible import log_infeasible_constraints
from decimal import Decimal
import matplotlib.pyplot as plt
pyomo.__version__
class Nonlinear_MPC:
    def __init__(self,x_lead_list,X_start,state_number,m,r,dt,H,a_min,a_max,f0,f1,f2,gamma_clf,gamma_cbf,slack_clf,slack_cbf,v_d,Th,cd,g,v0):
        
        # Model
        model = ConcreteModel()
        model.xk_num = RangeSet(1,H)
        model.uk_num = RangeSet(0,H)
        model.duk_num = RangeSet(0,H-1)
        
        # Parameters
        model.x_obj = x_lead_list
        model.dt = dt
        model.m = m
        model.state_number = state_number
        model.X0 = Param(RangeSet(0,state_number-1),initialize={0:X_start[0],1:X_start[1]},mutable = True)
        model.R = Param(RangeSet(1,5),initialize={1:10,2:10,3:1000,4:30},mutable=True)

        # Variables
        model.X = Var(RangeSet(0,state_number-1),model.uk_num)
        model.a = Var(model.uk_num,bounds=(a_min,a_max))
        model.slack_cbf = Var(model.uk_num,bounds=(0,1))
        
        # Constraint
        model.X0_update = Constraint(RangeSet(0,state_number-1),rule=lambda model,i:model.X[i,0]==model.X0[i])
        model.x_update = Constraint(model.uk_num,rule=lambda model,k:model.X[0,k+1] == model.X[0,k]+model.X[1,k]*model.dt
                                if k<=H-1 else Constraint.Skip)
        model.v_update = Constraint(model.uk_num,rule=lambda model,k:model.X[1,k+1] == model.X[1,k]+model.a[k]*model.dt/model.m-model.X[1,k]*(f1+2*f2*model.X[1,k])*model.dt/model.m
                                if k<=H-1 else Constraint.Skip)
#         model.clf_update = Constraint(model.uk_num,rule=lambda model,k:(v_d-model.X[1,k+1])-(v_d-model.X[1,k])+gamma_clf*(model.X[1,k]-v_d) <= model.slack_clf[k]
#                                 if k<=H-1 else Constraint.Skip)
        model.cbf_update = Constraint(model.uk_num,rule=lambda model,k:((model.x_obj[k+1]-model.X[0,k+1])**2-r**2)-((model.x_obj[k] - model.X[0,k])**2-r**2)+gamma_cbf*((model.x_obj[k]-model.X[0,k])**2-r**2)>= -model.slack_cbf[k]
                                if k<=H-1 else Constraint.Skip)
        
        #Objectives
        model.sbRsb = model.R[1]*sum(model.slack_cbf[i]**2 for i in model.uk_num)
        model.aRa = model.R[2]*sum(((model.a[i]-(f0+f1*model.X[1,i]+f2*model.X[1,i]**2))/model.m)**2 for i in model.uk_num)
        model.daRda = model.R[3]*sum((((model.a[i+1]-(f0+f1*model.X[1,i+1]+f2*model.X[1,i+1]**2))/model.m)-((model.a[i]-(f0+f1*model.X[1,i]+f2*model.X[1,i]**2))/model.m))**2 for i in model.duk_num)
        model.vRv = model.R[4]*sum(((model.X[1,i+1]-v_d))**2 for i in model.duk_num)
        model.obj = Objective(expr=model.sbRsb+model.vRv+model.aRa+model.daRda,sense=minimize)

        self.iN = model
    
    def update_lead(self,x_lead_list):
        self.iN.x_obj = x_lead_list
        
    def Solve(self,X_state):
        self.iN.X0.reconstruct({0:X_state[0], 1: X_state[1]})
        self.iN.x_update.reconstruct()
        self.iN.v_update.reconstruct()
#         self.iN.clf_update.reconstruct()
        self.iN.cbf_update.reconstruct()
        
        SolverFactory("ipopt").solve(self.iN)
        log_infeasible_constraints(self.iN)
        x_opt = [self.iN.X[0,k]() for k in self.iN.xk_num]
        v_opt = [self.iN.X[1,k]() for k in self.iN.xk_num]
        ax_opt = [self.iN.a[k]() for k in self.iN.uk_num]
        slack_cbf_opt = [self.iN.slack_cbf[k]() for k in self.iN.uk_num]
        return x_opt, v_opt, ax_opt, slack_cbf_opt
    
    
