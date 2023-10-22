#%% 1. Import Libraries
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sympy import bell
from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
import gurobipy as gp
from gurobipy import *
from gurobipy import tupledict
import time
import cmath
import math
#%% 2. Arrival Data
#periods of constant arrival
lambda_vec=[3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525]
lambda_vec = [i / 36 for i in lambda_vec]
# lambda_vec=[0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.0955153385708942,0.0955153385708942,0.0955153385708942,0.0955153385708942,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.00502712308267864,0.0502712308267864,0.0502712308267864,0.0502712308267864,0.0502712308267864]


#Crisis ArrivalRate
weekly_lam=[0.0078125,0.0129774305555556,0.0142361111111111,0.0173177083333333,0.0197048611111111,0.0225694444444444,0.0250868055555556,0.0246527777777778,0.0262586805555556,0.0273003472222222,0.0264756944444444,0.0238715277777778,0.0202690972222222,0.0180555555555556,0.0168836805555556,0.0133246527777778,0.00677083333333333,0.000477430555555556
]

lambda_vec=[]
for i in weekly_lam:
    lambda_vec=lambda_vec+[i]*40
    #%%
    #Crisis Arrival Rate
weekly_lam_c=[0.01042,0.0125,0.01944,0.01528,0.02083,0.02639,0.02292,0.02292,0.02847,0.01736,0.02431,0.02222,0.01597,0.00417,0.00208,0.00069,0.00069,0.00069
]
lambda_c=[]
for i in weekly_lam_c:
    lambda_c=lambda_c+[i]*40
    
M=100000
# Regular Arrival Rate
lambda_r=[3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,3.26025,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,4.181625,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.54375,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.217725,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.43035,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.161025,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,3.302775,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.877525,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.735775,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,2.35305,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.978075,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.269325,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.6237,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.212625,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.014175,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525,0.042525]
lambda_r = [i / 36 for i in lambda_r]
#%% 3. Debugging Block
lambda_vec=lambda_vec[1:40]
print(len(lambda_vec))
# len(lambda_vec)
#%% 4. Linearized Constraints: min prop
#Linearized
def minprop(lambda_vec):

      # Create the model within the Gurobi environment
      model = gp.Model("Test")
      
      s = model.addMVar(len(lambda_vec), vtype= GRB.BINARY,name="s")
      h = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="h")
      z = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="z")
      
      #OBJECTIVE
      model.setObjective(sum(s), GRB.MINIMIZE)



      #BOUNDARY CONDTIONS ON X AND H
      model.addConstr(h[-1]==1)
      model.addConstr(s[-1]==1)
      
         
      
      #NON LINEAR CONSTRAINTS
      # model.addConstrs(h[i]==s[i]+(1-s[i])*(1+h[i+1]) for i in range(len(lambda_vec)-1))
      
      
    # #LINEARIZING CONSTRAINTS
      model.addConstrs(h[i]==h[i+1]+1-z[i] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]<=s[i]*M for i in range(len(lambda_vec)))
      model.addConstrs(z[i]<=h[i+1] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]>=h[i+1]-M*(1-s[i]) for i in range(len(lambda_vec)-1))


      #UNDERUTILIZATION CONSTRAINT
      model.addConstrs(h[i]<=1/lambda_vec[i] for i in range(len(lambda_vec)))



      #SOLVE
      model.update()
      model.optimize()
      
      
      #RETURN AND PRINT
      print (h.x)
      print (s.x)
      return (sum(s.x))
      
      
      
#%%  5. MIN WAIT TIME
#Linearized
def minWT(lambda_vec,min_slots):
      M=999999

      # Create the model within the Gurobi environment
      model = gp.Model("Test")
      model.params.NonConvex = 2
      
      s = model.addMVar(len(lambda_vec), vtype= GRB.BINARY,name="s")
      h = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="h")
      z = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="z")
      zz =model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="zz")

      model.setObjective(sum(lambda_vec[i]*h[i]*zz[i] for i in range(len(lambda_vec))), GRB.MINIMIZE)




      #BOUNDARY CONDTIONS ON X AND H
      model.addConstr(h[-1]==1)
      model.addConstr(s[-1]==1)

      #NON LINEAR CONSTRAINTS
      # model.addConstrs(h[i]==s[i]+(1-s[i])*(1+h[i+1]) for i in range(len(lambda_vec)-1))
      
      
      # #LINEARIZING CONSTRAINTS
      model.addConstrs(h[i]==h[i+1]+1-z[i] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]<=s[i]*M for i in range(len(lambda_vec)))
      model.addConstrs(z[i]<=h[i+1] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]>=h[i+1]-M*(1-s[i]) for i in range(len(lambda_vec)-1))


      #UNDERUTILIZATION CONSTRAINT
      model.addConstrs(h[i]<=1/lambda_vec[i] for i in range(len(lambda_vec)))

      #MINPROP CONSTRAINT
      model.addConstr(sum(s)<=min_slots)
      
      #REMOVE DIVIDE TERM
      model.addConstrs(zz[i]*(1-lambda_vec[i]*h[i])==1 for i in range(len(lambda_vec))) 

      #Solve
      model.update()
      model.optimize()

      #return
      print(sum(s.x))
      # print (sum(a)
      print (h.x)
      print (s.x)
#%% 5. MIN WAIT TIME_NEW
#Linearized
def minWT_new(lambda_vec,min_slots):
      M=9999999
      # Create the model within the Gurobi environment
      model = gp.Model("Test")
      model.params.NonConvex = 2
      model.params.MIPFocus = 3
      
      
      s = model.addMVar(len(lambda_vec), vtype= GRB.BINARY,name="s")
      h = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="h")
      z = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="z")
      zz = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="zz")

      model.setObjective(sum(zz[i] for i in range(len(lambda_vec))), GRB.MINIMIZE)




      #BOUNDARY CONDTIONS ON X AND H
      model.addConstr(h[-1]==1)
      model.addConstr(s[-1]==1)

      #NON LINEAR CONSTRAINTS
      # model.addConstrs(h[i]==s[i]+(1-s[i])*(1+h[i+1]) for i in range(len(lambda_vec)-1))
      
      
      # #LINEARIZING CONSTRAINTS
      model.addConstrs(h[i]==h[i+1]+1-z[i] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]<=s[i]*M for i in range(len(lambda_vec)))
      model.addConstrs(z[i]<=h[i+1] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]>=h[i+1]-M*(1-s[i]) for i in range(len(lambda_vec)-1))


      #UNDERUTILIZATION CONSTRAINT
      model.addConstrs(h[i]<=1/lambda_vec[i] for i in range(len(lambda_vec)))

      #MINPROP CONSTRAINT
      model.addConstr(sum(s)<=min_slots)
      
      #REMOVE DIVIDE TERM
      model.addConstrs(zz[i]*2*(1-(lambda_vec[i]*h[i]))>=(lambda_vec[i]**2)*(h[i]**2) for i in range(len(lambda_vec))) 

      #Solve
      model.update()
      model.optimize()

      #return
      print(sum(s.x))
      # print (sum(a)
      print (h.x)
      print (s.x)          
#%%
if __name__ == '__main__':
    #start time
    start_time = time.time()
    minWT_new(lambda_vec,10)  
    #time it took
    print("--- %s mins ---" % ((time.time() - start_time)/60))

#%%
# minprop(lambda_c)   
