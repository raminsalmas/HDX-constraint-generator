#!/usr/bin/env python3

""" This script belongs to Borysik group in the Department of Chemistry of King's College London.
@author: Ramin Ekhteiari Salmas
"""

import numpy as np
import re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


#-------------------------------------------------------------------------------------
RFU_EXP = "Alpha_lac.txt" #import the exp. RFU file name

RFU_all_model = "RFU_values_all_Alpha_1.txt" #import the all_modelled_RFU file name

Kint = "K_int_Alpha.txt" # Kint file

#----------------------------------------------------------------------------------


with open( RFU_all_model, 'r') as rf2:
        fil = rf2.readlines()
        global RUN
        
        RUN=[ii.startswith("Rmse")for ii in fil].count(True)
        print("RUN:", RUN)



K_int = []

with open(Kint, 'r') as kin:
   for i in kin:
       if i.strip():
          kkk = i.split()
          K_int.append(float(kkk[2]))      
          #res_no.append(int(kkk[1]))
          #res_name.append(str(kkk[0]))
          #if kkk[0] == 'P':
             # Proline.append(int(kkk[1]))






with open(RFU_EXP, 'r') as read:
    for i, line in enumerate(read):
       if i == 0 and line.strip():
           global T
           T = [float(s) for s in re.findall(r'-?\d+\.?\d*', line)]

#Deterum_t = [[] for dtt in range(len(T))]
#RFU_exp_t = [[] for dtt in range(len(T))]
Deterum = [[] for dtt in range(RUN)]
RFU_exp = [[] for dtt in range(RUN)]

res_range = []
S_E_residue = []
peptide_No = []

#rfu= [[] for linei in range(50)]

bs= -1
with open( RFU_all_model, 'r') as rf2:
        for i, line in enumerate(rf2):
            sp = line.split()
            if len(sp)==1:
                bs = bs+1
                
            #if len(rmse_2) ==top_run +1 and len(sp)!=1:
            if len(sp)!=1:                             
               RFU_exp[bs].append([float(jli) for jli in sp]) 
              


np.random.seed(1980)
np.random.shuffle(RFU_exp)



with open(RFU_EXP, 'r') as read:
   
      
       for i, line in enumerate(read):                
           if i !=0 and line.strip():
               sp = line.split()
               exch = int(sp[-1]) - int(sp[-2])
               #print(exch)
               S_E_residue.append(sp[-2::])
               count = str(sp[0][1:]).count("P")
               if count !=0:
                   exch = exch - count
               peptide_No.append(int(sp[1]))
               res_range.append([int(Aoo) for Aoo in sp[-2:]])
               #print("ew", exch)
               #print(sp[0][1:])
               for jonas in range(RUN):
                  
                  RFU = RFU_exp[jonas][i-1]
                 # for ghh in range(len(T)):
                      #Deterum_t[ghh].append(float(sp[ghh+2])*exch)
                      #RFU_exp_t[ghh].append(float(sp[ghh+2]))                
                  data = [float(b)*exch for b in RFU]            
                  Deterum[jonas].append([g for g in data])
                  
                 
Amino_acid = []
Pep1 =[]
Pep2 =[]
for ii in range(len(res_range)):
    for jj in range(len(res_range) - ii):
        AB = np.absolute(np.subtract(res_range[ii], res_range[jj +ii]))
        #AC = np.absolute(np.subtract(AB[0],AB[1]))
        if AB[0]<=1 and AB[1]<=1 and AB[0] != AB[1]:
           #print(res_range[ii], res_range[jj +ii])
           #print(np.absolute(np.subtract(res_range[ii], res_range[jj +ii])))
           #print(ii+1, jj +ii+1)
           Pep1.append(ii)
           Pep2.append(jj +ii)
           if res_range[ii][0] != res_range[jj +ii][0]:
               #print(max(res_range[ii][0], res_range[jj +ii][0]))
               Amino_acid.append(max(res_range[ii][0], res_range[jj +ii][0]))
           else:
               #print(max(res_range[ii][1], res_range[jj +ii][1]))
               Amino_acid.append(max(res_range[ii][1], res_range[jj +ii][1]))
               

        #print(res_range[ii], res_range[jj +ii])
        #print(AB)
        #print(AC)
        
#print(Amino_acid)

               
Kobs = [[] for lsiss in range(RUN)]

pep_sele = []
for ijl in range(len(Pep1)):
    pep_sele.append(str([peptide_No[Pep1[ijl]],peptide_No[Pep2[ijl]]]))
#print(pep_sele)
#print(Pep1)

for onls in range(RUN):
   
   for lson in range(len(Amino_acid)):
       
       x= np.array(T)
       y = np.absolute(np.subtract(Deterum[onls][Pep1[lson]], Deterum[onls][Pep2[lson]]))
       

       def func(x, a):
           return 1-np.exp(-a*x)


       popt, pcov = curve_fit(func, x, y)

       Kobs[onls].append(popt[0])
       print(pep_sele[lson],Amino_acid[lson], "---", popt[0]) #screen
    
    #plt.plot(x, func(x, *popt))


res_rep = []
with open("output_RUN_%s" %(RFU_EXP), 'w') as out:
   out.write("Res. Peptides\n")
   for i, (Ami, pepp) in enumerate(zip(Amino_acid, pep_sele)):
         kobs_all = []
         
         if Ami not in res_rep:
                     if i !=0:
                        out.write("\n")

                     res_rep.append(Ami)

                     #out.write("%d %s  " %(Ami, pepp.replace(", ", '-')))
                     out.write("%d  " %(Ami))
                     for ll, kob in enumerate(Kobs):
                        #print(kob)
                        kobs_all.append(kob[i])
                        #out.write("%.15f " %(kob[i]))

                     #out.write("  %.15f" %(np.median(kobs_all)))
                     fi = np.log(np.divide(K_int[Ami],kob[i]))
                     fi1 = np.divide(K_int[Ami], np.exp(fi+2))
                     fi__1 = np.divide(K_int[Ami], np.exp(fi-2))
                     out.write("  %.15f %.15f" %(fi1, fi__1))
     





print(res_rep)



    
    
