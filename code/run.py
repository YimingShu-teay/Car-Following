import numpy as np
from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import math
from math import factorial
from pyomo.util.infeasible import log_infeasible_constraints
from decimal import Decimal
import matplotlib.pyplot as plt
from dynamical import *
from solver import Nonlinear_MPC
from IDM import IDM
from utils import *

#参数常量
Params = {}
Params['Nx'] = 2
Params['Nu'] = 1
Params['H'] = 5
Params['dt'] = 0.1
Params['m'] = 1650
Params['g'] = 9.81
Params['ca'] = 0.2
Params['cd'] = 0.2
Params['f0'] = 0.1
Params['f1'] = 1.8
Params['f2'] = 0.15
Params['Line_Y'] = 10.0
Params['gamma_clf'] = 0.5
Params['gamma_cbf'] = 0.2
Params['slack_clf'] = 2e-2
Params['slack_cbf'] = 2e-2
Params['a_max'] = Params['ca']*Params['m']*Params['g']
Params['a_min'] = -Params['cd']*Params['m']*Params['g']
Params['v_d'] = 25.0
Params['x_lead0'] = 30.0
Params['y_lead0'] = 10.0
Params['v_lead0'] = 15.0
Params['a_lead0'] = 0.0
Params['A_max'] = 1.0
Params['V0'] = 15.0
Params['X_start'] = [5.0,10.0]
Params['MAX_TIME'] = 10.0
Params['r'] = 10.0
Params['Line_Y'] = 10.0
Params['Th'] = 1.8
Params['v0'] = 14.0
Params['show_animation'] = True


def run_simulation(x_lead0,y_lead0,v_lead0,a_lead0,X_start,H,dt,MAX_TIME,show_animation,Nx,r,m,a_min,a_max,f0,f1,f2,
                   gamma_clf,gamma_cbf,slack_clf,slack_cbf,v_d,A_max,V0,Th,cd,g,v0):
    #list用来记录，之后用来画图
    time = 0.0
    t = [0.0]
    
    ugv = KinematicModel(X_start[0],X_start[1],dt,f0,f1,f2,m)
    ugv_lead = IDM(v_lead0,x_lead0,y_lead0,a_lead0,A_max,dt,V0)
    x_lead_list = ugv_lead.get_x_lead_list(H)
    
    x = [ugv.x]
    v = [ugv.v]
    ax = [0.0]
    x_lead = [ugv_lead.Leader_X]
    v_lead = [ugv_lead.Leader_V]
    a_lead = [ugv_lead.Leader_a]
    slack = [0.0]
    
    while MAX_TIME >= time:
        X0 = [ugv.x, ugv.v]
        x_lead_list = ugv_lead.get_x_lead_list(H)
        MPC_Model = Nonlinear_MPC(x_lead_list,X0,Nx,m,r,dt,H,a_min,a_max,f0,f1,f2,gamma_clf,gamma_cbf,slack_clf,slack_cbf,v_d,Th,cd,g,v0)
        MPC_Model.update_lead(x_lead_list)
        x_opt, v_opt, ax_opt, slack_clf_opt = MPC_Model.Solve(X0)
        slack.append(slack_clf_opt)
        if ax_opt is not None:
            ai = ax_opt[0]
 
        ugv.update_state(ai)
        ugv_lead.update_state()
        print(ugv_lead.Leader_X,ugv.x,Th,ugv.v)
        print(ugv_lead.Leader_X - ugv.x-10.0)
        time = time + dt
        t.append(time)
        x.append(ugv.x)
        v.append(ugv.v)
        ax.append(ai)
        
        x_lead.append(ugv_lead.Leader_X)
        v_lead.append(ugv_lead.Leader_V)
        a_lead.append(ugv_lead.Leader_a)
        
        if show_animation:
            
            plt.cla()
            plot_car(ugv.x, 10.0, 0, 0)
            plot_car(ugv_lead.Leader_X, ugv_lead.Leader_Y, 0, 0)
         
            plt.axis("equal")
            plt.grid(True)
            
            plt.pause(0.0001)

    return x,v,ax,x_lead,v_lead,a_lead,t

if __name__ == '__main__':
    x,v,ax,x_lead,v_lead,a_lead,t = run_simulation(Params['x_lead0'],Params['y_lead0'],Params['v_lead0'],Params['a_lead0'],
                   Params['X_start'],Params['H'],Params['dt'],Params['MAX_TIME'],Params['show_animation'],
                   Params['Nx'],Params['r'],Params['m'],Params['a_min'],Params['a_max'],
                   Params['f0'],Params['f1'],Params['f2'],Params['gamma_clf'],Params['gamma_cbf'],
                   Params['slack_clf'],Params['slack_cbf'],Params['v_d'],
                   Params['A_max'],Params['V0'],Params['Th'],Params['cd'],Params['g'],Params['v0'])
    plt.subplots()
    plt.plot(t, v, "-r", label="speed")
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.ylabel("Speed [m/s]")

    plt.subplots()
    plt.plot(t, ax, "-r", label="speed")
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.ylabel("Acc [m/s]")
    
    plt.show()