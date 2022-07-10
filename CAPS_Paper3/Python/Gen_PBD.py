#%% 
# 1. IMPORT LIBRARIES
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sympy import bell
from tqdm import tqdm
import numpy as np 
# import scipy.stats as stats
# import sys
import pandas as pd
from gurobipy import *
from gurobipy import tupledict
# import time
import cmath
import math
#%%
# 2. POISSON BINOMIAL PDF
pi=math.pi
compuv=complex(0,1)

risk_vec=[0.00501450824706418,0.00501450824706418,0.00501450824706418,0.00501450824706418,0.139992090324993,0.139992090324993,0.139992090324993,0.139992090324993,0.0910955797481352,0.0910955797481352,0.0910955797481352,0.0910955797481352,0.178035276440403,0.178035276440403,0.178035276440403,0.178035276440403,0.0490285432565714,0.0490285432565714,0.0490285432565714,0.0490285432565714]
blks=len(risk_vec)
def PBD(Node_1,Node_2,m):
    pdf=0.0
    cdf=0.0
    for n in range(0,Node_2-Node_1+1):
        prod=1.0
        for k in range(Node_1,Node_2):
            prod=prod*(risk_vec[k]*cmath.exp(compuv*2*pi*n/(Node_2-Node_1+1.0))+(1-risk_vec[k]))
        pdf=pdf+cmath.exp(-compuv*2*pi*n*m/(Node_2-Node_1+1.0))*prod 
        if n==0:
            cdf=cdf
        else:
            cdf=cdf+(1.0-cmath.exp(-compuv*2*pi*n*m/(Node_2-Node_1+1.0)))/(1.0-cmath.exp(-compuv*2*pi*n/(Node_2-Node_1+1.0)))*prod           
    return (pdf/(Node_2-Node_1+1.0)).real, (1.0-float(m)/(Node_2-Node_1+1.0)-cdf/(Node_2-Node_1+1.0)).real
#%% create P_vec
P_vec= np.zeros((blks+1,blks))
P_vec= [[PBD(0,i,j)[0] for i in range(blks)] for j in range(blks+1)] 
# %%
