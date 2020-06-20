#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 22:24:23 2020

@author: jguterl
"""
import sys,os
#import omfit 
#import ImportEqFiles
Shot=179846
#Out=ImportEqFiles.ImportEqFiles(167196,3000,ForceReload=False)


os.chdir( '/Users/torvaltz/Desktop/INGRID/src')
import Ingrid
from Ingrid import *
InputFile='/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS1_modif.yml'
WTargetFile='/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS_odt.txt'
ETargetFile='/Users/torvaltz/Desktop/INGRID/data/SNL/USN/itp4.txt'
EqFile='/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/g175816.03000'
igrd=Ingrid(yaml=InputFile,w1=WTargetFile,e1=ETargetFile,eq=EqFile)

igrd.yaml['grid_params']['grid_generation']['radial_f_core']='x, x**1.4'
igrd.yaml['grid_params']['grid_generation']['radial_f_sol']='x,(1-exp(-(x)/0.50))/(1-exp(-1/0.50))'
igrd.yaml['grid_params']['grid_generation']['radial_f_pf']='x,1-(1-exp(-(1-x)/0.90))/(1-exp(-1/0.90))'
igrd.yaml['target_plates']['plate_W1']['poloidal_f']='x, 1-(1-exp(-(1-x)/0.30))/(1-exp(-1/0.30))'
igrd.yaml['target_plates']['plate_E1']['poloidal_f']='x, (1-exp(-(x)/0.30))/(1-exp(-1/0.30))'

igrd.yaml['target_plates']['plate_W1']['np_local']=4
igrd.yaml['target_plates']['plate_E1']['np_local']=4

igrd.yaml['grid_params']['grid_generation']['np_sol']=4
igrd.yaml['grid_params']['grid_generation']['np_core']=4
igrd.yaml['grid_params']['grid_generation']['np_pf']=4

igrd.yaml['grid_params']['grid_generation']['nr_sol']=4
igrd.yaml['grid_params']['grid_generation']['nr_core']=4
igrd.yaml['grid_params']['grid_generation']['nr_pf']=4


igrd.Setup()
igrd.ShowSetup()
#igrd.current_topology.construct_patches()
#igrd.current_topology.SavePatches('/home/jguterl/Dropbox/python/GridGenerator/D3DSAS/patches.yml')

igrd.current_topology.LoadPatches('/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/patches.yml')
igrd.current_topology.patch_diagram()
igrd.current_topology.CheckPatches()
self=igrd.current_topology.patches

#fig=plt.figure()
#ax=plt.gca()
plt.show()
plt.draw()
plt.ion()
# List=self.patches
# for patch in List:
#             (nr_cells,np_cells)=self.GetNpoints(patch)
#             (_radial_f,_poloidal_f)=self.GetFunctions(patch)
#             print('# Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(patch.patchName,np_cells,nr_cells,inspect.getsource(_poloidal_f),inspect.getsource(_radial_f)))
#             patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f)
#             self.AdjustPatch(patch)
#             patch.plot_subgrid()
#             plt.pause(1)

#             #input("Press Enter to continue...")
#             List.remove(patch)

pf1={'poloidal_f':'x,(x)**1.5','np_cells':2}
pf2={'poloidal_f':'x,1-(1-x)**1.5','np_cells':2}
#pf1={'poloidal_f':'x,(x)**1','np_cells':10}
#pf2={'poloidal_f':'x,1-(1-x)**1','np_cells':10}
Dic={'OCB':pf1,'OSB':pf1,'ICB':pf2,'ISB':pf2,'IST':{'np_cells':2},'ICT':{'np_cells':2},'OST':{'np_cells':2},'OCT':{'np_cells':2}}
for p in igrd.current_topology.patches.values():
    p.Verbose=True
CorrectionDistortion={'ThetaMin':70,'ThetaMax':110,'Resolution':1000,'Active':False}
igrd.current_topology.CorrectDistortion={'OSB':CorrectionDistortion,'OCB':CorrectionDistortion,'all':CorrectionDistortion}
igrd.current_topology.Verbose=True
plt.ion()
plt.show()
igrd.Verbose=True
igrd.current_topology.patches['IDL'].Verbose=True
igrd.current_topology.patches['IPF'].Verbose=True
igrd.CreateSubgrid(ShowVertices='detailed',RestartScratch=True,NewFig=False,OptionTrace='theta',ExtraSettings=Dic,ListPatches='all')
['ISB','ICB']
#%%
#igrd.CreateSubgrid(ShowVertices=True,RestartScratch=False,NewFig=False)
igrd.current_topology.gridue_params
fig,ax=plt.subplots()
Plotgridue(igrd.current_topology.gridue_params,ax=plt.gca(),Verbose=True,facecolor=None,edgecolor='blue')
FN='/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/'+'gridue_d3dSAS_test'
igrd.export(fname=FN)
M=UEDGEMesh()
g=igrd.current_topology.gridue_params
M.SetGrid(g)
M.ShowMesh(ax)
ListCell=CheckOverlapCells(igrd.current_topology.gridue_params)
for (p1,p2) in ListCell:
    M.ShowCell(p1)
    M.ShowCell(p2)
#igrd.StartGUI()

x=np.linspace(0,1,20)
f=lambda x:(1-(1-x)**2)-0.005*(exp(-x*10)-1)*(exp(-(1-x)*10000)-1)
fig,ax=plt.subplots()
ax.plot(x,f(x),marker='o')
ax.plot(x,(1-(1-x)**2),marker='o')

x=np.linspace(0,1,20)
a=0.10
f=lambda x:(1-exp(-(x)/a))/(1-exp(-1/a))
fig,ax=plt.subplots()
#ax.plot(x,f(x),marker='o')
ax.plot(x,(1-exp(-(x)/0.30))/(1-exp(-1/0.30)),marker='o')
ax.plot(x,1-(1-exp(-(1-x)/0.30))/(1-exp(-1/0.30)),marker='s')