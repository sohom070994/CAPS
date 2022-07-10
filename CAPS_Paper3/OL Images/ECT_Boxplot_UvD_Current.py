# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%%

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
 
#%%
path="C:\\Users\\Sohom\\Desktop\\CAPS_Paper3\\OL Images\\Matplotlib_Images\\"
Dyn_8_4=np.array([1.29454015100788,1.06397323914981,0.850406732899272,0.978470800571103,1.66638519989723,1.14221851220086,0.884001418058552,0.705241705915138,1.47343476206528,0.946395750020106,1.02957775280834,1.76041358175043,0.716787404089961,0.952999941293177,0.786460689421804,0.806047448666023,0.971088358652947,1.32183421550266,1.08327241664595,1.06065377637775,1.24933606698728,0.92684665902443,1.12904772443327,0.895348654462842,0.805308328481748,1.40681187619992,0.936669339399125,1.5751405933754,1.10207601498021,1.01007019755251,1.03167765966828,0.83408236698228])
Uni_8_4=np.array([1.06330566996769,1.87203159774816,1.82225232737185,2.10250267913331,3.07395093011277,1.09142783659274,1.11524950911018,1.89063766623075,2.91315606029831,1.35869151951855,2.05314277971196,1.94359777261861,0.845843047253362,1.19873580673075,1.63031303979733,1.25782373809388,2.32316302439353,1.25965102768721,3.38204005663005,2.68848663930949,2.33073349687167,0.953798485237577,0.937083401946218,1.45551965830843,0.925307891589898,1.8138364979044,1.39921219464888,1.3344236638953,1.81328380164214,1.37352908644395,1.52866418243955,0.717983286635782])
#%%
data=[Uni_8_4,Dyn_8_4]

fig = plt.figure(figsize=(6,4))

#change fatness
xl=0  
xu=0.75

#change distance from ends
s=3.5

# Creating axes instance
ax = fig.add_axes([0.1, 0.1, .85, .85])#if stating two args 1 then axis will not show

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
# Creating plot
bp = plt.boxplot(data,patch_artist=True,positions=[(xu-xl)/s,xu-(xu-xl)/s],showmeans=True,showfliers=False,
                 meanprops={"marker":"x","markeredgecolor":"black"},
                 medianprops=dict(linestyle='dashed', linewidth=1, color='black')
                 )


colors = ['red', 'green']
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)


plt.xlim(xl,xu)
plt.grid(axis = 'y',alpha=0.3)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) 

# plt.xlabel('categories')
plt.ylabel('Expected crisis time (hrs)')

ax.legend(bp["boxes"], ["Uniform", "Dynamic"], loc='upper right')

fig1 = plt.gcf()

# show plot
# plt.plot()
#%%
fig1.savefig(path+'ECT-4pct_BP.pdf',dpi=100)  
fig1.show()