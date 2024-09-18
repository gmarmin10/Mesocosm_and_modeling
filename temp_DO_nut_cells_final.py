'''
Created on Jun 15, 2022

average, st.dev for daily temperature, pH, and dissolved oxygen

@author: garmin
'''
from numpy import *
from matplotlib.pyplot import plot, show, scatter, legend, xlim, ylim, figure, xlabel, xticks, ylabel, yticks,errorbar, title, ticklabel_format


#eukaryotic tanks
data=genfromtxt('avg_temp_pH_DO.csv',delimiter=',') 
day = data[2:,0]

tank1 = data[2:,1:13]
T_1 = tank1[:,0]
TU_1 = tank1[:,2]
TL_1 = tank1[:,3]
error_T_1 = [(T_1-TL_1),(-T_1+TU_1)]
DO_1 = tank1[:,8]
DOU_1 = tank1[:,10]
DOL_1 = tank1[:,11]
error_DO_1 = [(DO_1-DOL_1),-(DO_1-DOU_1)]


tank5 = data[2:,13:25]
T_5 = tank5[:,0]
TU_5 = tank5[:,2]
TL_5 = tank5[:,3]
error_T_5 = [(T_5-TL_5),(-T_5+TU_5)]
DO_5 = tank5[:,8]
DOU_5 = tank5[:,10]
DOL_5 = tank5[:,11]
error_DO_5 = [(DO_5-DOL_5),-(DO_5-DOU_5)]


tank9 = data[2:,25:37]
T_9 = tank9[:,0]
TU_9 = tank9[:,2]
TL_9 = tank9[:,3]
error_T_9 = [(T_9-TL_9),(-T_9+TU_9)]
DO_9 = tank9[:,8]
DOU_9 = tank9[:,10]
DOL_9 = tank9[:,11]
error_DO_9 = [(DO_9-DOL_9),-(DO_9-DOU_9)]


tank2 = data[2:,37:49]
T_2 = tank2[:,0]
TU_2 = tank2[:,2]
TL_2 = tank2[:,3]
error_T_2 = [(T_2-TL_2),(-T_2+TU_2)]
DO_2 = tank2[:,8]
DOU_2 = tank2[:,10]
DOL_2 = tank2[:,11]
error_DO_2 = [(DO_2-DOL_2),-(DO_2-DOU_2)]


tank6 = data[2:,49:61]
T_6 = tank6[:,0]
TU_6 = tank6[:,2]
TL_6 = tank6[:,3]
error_T_6 = [(T_6-TL_6),(-T_6+TU_6)]
DO_6 = tank6[:,8]
DOU_6 = tank6[:,10]
DOL_6 = tank6[:,11]
error_DO_6 = [(DO_6-DOL_6),-(DO_6-DOU_6)]


tank12 = data[2:,61:73]
T_12 = tank12[:,0]
TU_12 = tank12[:,2]
TL_12 = tank12[:,3]
error_T_12 = [(T_12-TL_12),(-T_12+TU_12)]
DO_12 = tank12[:,8]
DOU_12 = tank12[:,10]
DOL_12 = tank12[:,11]
error_DO_12 = [(DO_12-DOL_12),-(DO_12-DOU_12)]


figure(1, figsize=(8,6.5))
scatter(day,T_1, color='#88CCEE', s=60, label='RE-1')
errorbar(day,T_1,yerr=error_T_1,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_5, color='#228833', s=60,  label='RE-5')
errorbar(day,T_5,yerr=error_T_5,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_9, color='#4477AA', s=60, label='RE-9')
errorbar(day,T_9,yerr=error_T_9,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_2, color='#DDCC77', s=60, marker='s', label='HE-2')
errorbar(day,T_2,yerr=error_T_2,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_6, color='#CC6677', s=60,marker='s', label='HE-6')
errorbar(day,T_6,yerr=error_T_6,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_12, color='#882255', s=60,marker='s', label= 'HE-12')
errorbar(day,T_12,yerr=error_T_12,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
title('Eukaryotic Tanks', fontsize=20)
#xlabel('Day of Experiment', fontsize=16)
#ticklabel_format(axis="y", style="sci", scilimits=(0,0))
xticks(fontsize=16)
ylabel('Water Temperature ($^o$C)', fontsize=18)
yticks(fontsize=16)
ylim(15,30)
#legend(fontsize='18',frameon=False,ncol=2) 
legend(loc='upper center', ncol=3,fontsize='large',frameon=False)       

 
figure(2, figsize=(8,6.5))
scatter(day,DO_1, color='#88CCEE', s=40, label='RE-1')
errorbar(day,DO_1,yerr=error_DO_1,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_5, color='#228833', s=40,  label='RE-5')
errorbar(day,DO_5,yerr=error_DO_5,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_9, color='#4477AA', s=40,  label='RE-9')
errorbar(day,DO_9,yerr=error_DO_9,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_2, color='#DDCC77', s=40, marker='s', label='HE-2')
errorbar(day,DO_2,yerr=error_DO_2,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_6, color='#CC6677', s=40, marker='s', label=' HE-6')
errorbar(day,DO_6,yerr=error_DO_6,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_12, color='#882255', s=40,  marker='s',label= 'HE-12')
errorbar(day,DO_12,yerr=error_DO_12,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
xticks(fontsize=12)
title('Eukaryotic Tanks', fontsize=20)
ylabel('DO (mg L$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(0,18)

data = genfromtxt('eukaryotic_nutrients.csv',delimiter=',')

day = data[1:,0]


no3_1 = data[1:,2]
TDN_1 = data[1:,3]
TN_1 = data[1:,4]
react_P_1 = data[1:,5]
TDP_1 = data[1:,6]
TP_1 = data[1:,7]

no3_5 = data[1:,9]
TDN_5 = data[1:,10]
TN_5 = data[1:,11]
react_P_5 = data[1:,12]
TDP_5 = data[1:,13]
TP_5 = data[1:,14]

no3_9 = data[1:,16]
TDN_9 = data[1:,17]
TN_9 = data[1:,18]
react_P_9 = data[1:,19]
TDP_9 = data[1:,20]
TP_9 = data[1:,21]

no3_2 = data[1:,23]
TDN_2 = data[1:,24]
TN_2 = data[1:,25]
react_P_2 = data[1:,26]
TDP_2 = data[1:,27]
TP_2 = data[1:,28]

no3_6 = data[1:,30]
TDN_6 = data[1:,31]
TN_6 = data[1:,32]
react_P_6 = data[1:,33]
TDP_6 = data[1:,34]
TP_6 = data[1:,35]

no3_12 = data[1:,37]
TDN_12 = data[1:,38]
TN_12 = data[1:,39]
react_P_12 = data[1:,40]
TDP_12 = data[1:,41]
TP_12 = data[1:,42]

TDN_1= TDN_1/(14*1e6)
TDN_5=TDN_5/(14*1e6)
TDN_9= TDN_9/(14*1e6)
TDN_2=TDN_2/(14*1e6)
TDN_6= TDN_6/(14*1e6)
TDN_12=TDN_12/(14*1e6)


######## Total Dissolved Nitrogen #######
figure(3, figsize=(8,6.5))
scatter(day,TDN_1*1e6, color='#88CCEE', s=40, label='RE-1')
scatter(day,TDN_5*1e6, color='#117733', s=40, label='RE-5')
scatter(day,TDN_9*1e6, color='#4B0082', s=40, label='RE-9')
scatter(day,TDN_2*1e6, color='#DDCC77',marker='s', s=40, label='HE-2')
scatter(day,TDN_6*1e6, color='#CC6677',marker='s', s=40, label=' HE-6')
scatter(day,TDN_12*1e6, color='#882255',marker='s', s=40, label= 'HE-12')
#xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
ylabel('TDN (\u03bcmol/L)', fontsize=16)
yticks(fontsize=12)
ylim(0,300)


TDP_1= TDP_1/(30.97*1e6)
TDP_5=TDP_5/(30.97*1e6)
TDP_9= TDP_9/(30.97*1e6)
TDP_2=TDP_2/(30.97*1e6)
TDP_6= TDP_6/(30.97*1e6)
TDP_12=TDP_12/(30.97*1e6)

##### Total Dissolved Phosphorus #######
figure(4, figsize=(8,6.5))
scatter(day,TDP_1*1e6, color='#88CCEE', s=40, label='RE-1')
scatter(day,TDP_5*1e6, color='#117733', s=40, label='RE-5')
scatter(day,TDP_9*1e6, color='#4B0082', s=40, label='RE-9')
scatter(day,TDP_2*1e6, color='#DDCC77',marker='s', s=40, label='HE-2')
scatter(day,TDP_6*1e6, color='#CC6677',marker='s', s=40, label=' HE-6')
scatter(day,TDP_12*1e6, color='#882255',marker='s', s=40, label= 'HE-12')
xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
ylabel('TDP (\u03bcmol/L)', fontsize=16)
yticks(fontsize=12)
ylim(0,3)

data=genfromtxt('Cell_density_final_eukaryotic.csv',delimiter=',') 

day = data[2:17,1]

### cell density ###

T1_av = data[2:17,3]
T1_upper = data[2:17,4]
T1_lower = data[2:17,5]
error_T1 = [(T1_av-T1_lower),-(T1_av-T1_upper)]

T2_av = data[2:17,7]
T2_upper = data[2:17,8]
T2_lower = data[2:17,9]
error_T2 = [(T2_av-T2_lower),-(T2_av-T2_upper)]

T5_av = data[2:17,10]
T5_upper = data[2:17,12]
T5_lower = data[2:17,13]
error_T5 = [(T5_av-T5_lower),-(T5_av-T5_upper)]

T6_av = data[2:17,14]
T6_upper = data[2:17,16]
T6_lower = data[2:17,17]
error_T6 = [(T6_av-T6_lower),-(T6_av-T6_upper)]

T9_av = data[2:17,18]
T9_upper = data[2:17,20]
T9_lower = data[2:17,21]
error_T9 = [(T9_av-T9_lower),-(T9_av-T9_upper)]

T12_av = data[2:17,22]
T12_upper = data[2:17,24]
T12_lower = data[2:17,25]
error_T12 = [(T12_av-T12_lower),-(T12_av-T12_upper)]

##### Growth rate data import ####

grow=genfromtxt('Growth_rate.csv',delimiter=',')

T1 = grow[1:16,1]
T2 = grow[1:16,2]
T5 = grow[1:16,3]
T6 = grow[1:16,4]
T9 = grow[1:16,5]
T12 = grow[1:16,6]


# convert cells/mL to cells/L
T1_av=T1_av*1000
T5_av=T5_av*1000
T9_av=T9_av*1000
T2_av=T2_av*1000
T6_av=T6_av*1000
T12_av=T12_av*1000

figure(5, figsize=(8,6.5))
scatter(day,T1_av, color='#88CCEE', s=60, label='RE-1')
errorbar(day,T1_av,yerr=error_T1,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T5_av, color='#228833', s=60,  label='RE-5')
errorbar(day,T5_av,yerr=error_T5,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T9_av, color='#4477AA', s=60,  label='RE-9')
errorbar(day,T9_av,yerr=error_T9,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T2_av, color='#DDCC77', s=60, marker='s', label='HE-2')
errorbar(day,T2_av,yerr=error_T2,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T6_av, color='#CC6677', s=60, marker='s', label=' HE-6')
errorbar(day,T6_av,yerr=error_T6,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T12_av, color='#882255', s=60,  marker='s',label= 'HE-12')
errorbar(day,T12_av,yerr=error_T12,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
xlabel('Day of Experiment', fontsize=18)
ticklabel_format(axis="y", style="sci", scilimits=(0,0))
xticks(fontsize=16)
ylabel('Cell density (cells L$^{-1}$)', fontsize=18)
yticks(fontsize=16)
ylim(top=5e8)



#mixed population tanks
data=genfromtxt('avg_temp_pH_DO.csv',delimiter=',') 

day = data[2:,0]


tank3 = data[2:,73:85]
T_3 = tank3[:,0]
TU_3 = tank3[:,2]
TL_3 = tank3[:,3]
error_T_3 = [(T_3-TL_3),-(T_3-TU_3)]
DO_3 = tank3[:,8]
DOU_3 = tank3[:,10]
DOL_3 = tank3[:,11]
error_DO_3 = [(DO_3-DOL_3),-(DO_3-DOU_3)]


tank8 = data[2:,85:97]
T_8 = tank8[:,0]
TU_8 = tank8[:,2]
TL_8 = tank8[:,3]
error_T_8 = [(T_8-TL_8),-(T_8-TU_8)]
DO_8 = tank8[:,8]
DOU_8 = tank8[:,10]
DOL_8 = tank8[:,11]
error_DO_8 = [(DO_8-DOL_8),-(DO_8-DOU_8)]


tank10 = data[2:,97:109]
T_10 = tank10[:,0]
TU_10 = tank10[:,2]
TL_10 = tank10[:,3]
error_T_10 = [(T_10-TL_10),-(T_10-TU_10)]
DO_10 = tank10[:,8]
DOU_10 = tank10[:,10]
DOL_10 = tank10[:,11]
error_DO_10 = [(DO_10-DOL_10),-(DO_10-DOU_10)]


tank4 = data[2:,109:121]
T_4 = tank4[:,0]
TU_4 = tank4[:,2]
TL_4 = tank4[:,3]
error_T_4 = [(T_4-TL_4),-(T_4-TU_4)]
DO_4 = tank4[:,8]
DOU_4 = tank4[:,10]
DOL_4 = tank4[:,11]
error_DO_4 = [(DO_4-DOL_4),-(DO_4-DOU_4)]


tank7 = data[2:,121:133]
T_7 = tank7[:,0]
TU_7 = tank7[:,2]
TL_7 = tank7[:,3]
error_T_7 = [(T_7-TL_7),-(T_7-TU_7)]
DO_7 = tank7[:,8]
DOU_7 = tank7[:,10]
DOL_7 = tank7[:,11]
error_DO_7 = [(DO_7-DOL_7),-(DO_7-DOU_7)]


tank11 = data[2:,133:145]
T_11 = tank11[:,0]
TU_11 = tank11[:,2]
TL_11 = tank11[:,3]
error_T_11 = [(T_11-TL_11),-(T_11-TU_11)]
DO_11 = tank11[:,8]
DOU_11 = tank11[:,10]
DOL_11 = tank11[:,11]
error_DO_11 = [(DO_11-DOL_11),-(DO_11-DOU_11)]

figure(6, figsize=(8,6.5))
scatter(day,T_3, color='#88CCEE', s=60, marker='P' ,label='RM-3')
errorbar(day,T_3,yerr=error_T_3,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_8, color='#228833', s=60, marker='P', label='RM-8')
errorbar(day,T_8,yerr=error_T_8,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_10, color='#4477AA', s=60, marker='P', label='RM-10')
errorbar(day,T_10,yerr=error_T_10,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_4, color='#DDCC77', s=60, marker='v', label='HM-4')
errorbar(day,T_4,yerr=error_T_4,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_7, color='#CC6677', s=60, marker='v', label='HM-7')
errorbar(day,T_7,yerr=error_T_7,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T_11, color='#882255', s=60, marker='v', label= 'HM-11')
errorbar(day,T_11,yerr=error_T_11,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
title('Mixed Population Tanks', fontsize=20)
xticks(fontsize=16)
yticks(fontsize=16)
legend(loc='upper center', ncol=3,fontsize='large',frameon=False)       
ylim(15,30)

figure(7, figsize=(8,6.5))
scatter(day,DO_3, color='#88CCEE', s=40, marker='P', label='RM-3')
errorbar(day,DO_3,yerr=error_DO_3,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_8, color='#228833', s=40, marker='P', label='RM-8')
errorbar(day,DO_8,yerr=error_DO_8,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_10, color='#4477AA', s=40, marker='P', label='RM-10')
errorbar(day,DO_10,yerr=error_DO_10,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_4, color='#DDCC77', s=40, marker='v', label='HM-4')
errorbar(day,DO_4,yerr=error_DO_4,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_7, color='#CC6677', s=40, marker='v', label=' HM-7')
errorbar(day,DO_7,yerr=error_DO_7,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,DO_11, color='#882255', s=40, marker='v', label= 'HM-11')
errorbar(day,DO_11,yerr=error_DO_11,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
xticks(fontsize=12)
yticks(fontsize=12)
title('Mixed Population Tanks', fontsize=20)
ylim(0,18)


data = genfromtxt('prokaryotic_nutrients.csv',delimiter=',')

day = data[1:,0]

no3_3 = data[1:,2]
TDN_3 = data[1:,3]
TN_3 = data[1:,4]
react_P_3 = data[1:,5]
TDP_3 = data[1:,6]
TP_3 = data[1:,7]

no3_8 = data[1:,9]
TDN_8 = data[1:,10]
TN_8 = data[1:,11]
react_P_8 = data[1:,12]
TDP_8 = data[1:,13]
TP_8 = data[1:,14]

no3_10 = data[1:,16]
TDN_10 = data[1:,17]
TN_10 = data[1:,18]
react_P_10 = data[1:,19]
TDP_10 = data[1:,20]
TP_10 = data[1:,21]

no3_4 = data[1:,23]
TDN_4 = data[1:,24]
TN_4 = data[1:,25]
react_P_4 = data[1:,26]
TDP_4 = data[1:,27]
TP_4 = data[1:,28]

no3_7 = data[1:,30]
TDN_7 = data[1:,31]
TN_7 = data[1:,32]
react_P_7 = data[1:,33]
TDP_7 = data[1:,34]
TP_7 = data[1:,35]

no3_11 = data[1:,37]
TDN_11 = data[1:,38]
TN_11 = data[1:,39]
react_P_11 = data[1:,40]
TDP_11 = data[1:,41]
TP_11 = data[1:,42]

TDN_3=TDN_3/(14*1e6)
TDN_8=TDN_8/(14*1e6)
TDN_10=TDN_10/(14*1e6)
TDN_4=TDN_4/(14*1e6)
TDN_7=TDN_7/(14*1e6)
TDN_11=TDN_11/(14*1e6)

TDP_3=TDP_3/(30.97*1e6)
TDP_8=TDP_8/(30.97*1e6)
TDP_10=TDP_10/(30.97*1e6)
TDP_4=TDP_4/(30.97*1e6)
TDP_7=TDP_7/(30.97*1e6)
TDP_11=TDP_11/(30.97*1e6)


######## Total Dissolved Nitrogen #######
figure(8, figsize=(8,6.5))
scatter(day,TDN_3*1e6, color='#88CCEE',marker='P', s=40, label='RM-3')
scatter(day,TDN_8*1e6, color='#117733', marker='P',s=40, label='RM-8')
scatter(day,TDN_10*1e6, color='#4B0082', marker='P',s=40, label='RM-10')
scatter(day,TDN_4*1e6, color='#DDCC77', s=40,marker='v', label='HM-4')
scatter(day,TDN_7*1e6, color='#CC6677', s=40, marker='v',label=' HM-7')
scatter(day,TDN_11*1e6, color='#882255', s=40,marker='v', label= 'HM-11')
xticks(fontsize=12)
ylabel('TDN (\u03bcmol/L)', fontsize=16)
yticks(fontsize=12)
ylim(0,300)
 


##### Total Dissolved Phosphorus #######
figure(9, figsize=(8,6.5))
scatter(day,TDP_3*1e6, color='#88CCEE',marker='P', s=40, label='RM-3')
scatter(day,TDP_8*1e6, color='#117733',marker='P', s=40, label='RM-8')
scatter(day,TDP_10*1e6, color='#4B0082',marker='P', s=40, label='RM-10')
scatter(day,TDP_4*1e6, color='#DDCC77', s=40,marker='v', label='HM-4')
scatter(day,TDP_7*1e6, color='#CC6677', s=40, marker='v',label='HM-7')
scatter(day,TDP_11*1e6, color='#882255', s=40, marker='v',label= 'HM-11')
xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
ylabel('TDP (\u03bcmol/L)', fontsize=16)
yticks(fontsize=12)
ylim(0,3)



data=genfromtxt('cell_density_final_contaminant.csv',delimiter=',') 

day = data[2:17,1]

### cell density ###

T3_av = data[2:17,3]
T3_upper = data[2:17,4]
T3_lower = data[2:17,5]
error_T3 = [(T3_av-T3_lower),-(T3_av-T3_upper)]

T4_av = data[2:17,7]
T4_upper = data[2:17,8]
T4_lower = data[2:17,9]
error_T4 = [(T4_av-T4_lower),-(T4_av-T4_upper)]

T8_av = data[2:17,11]
T8_upper = data[2:17,12]
T8_lower = data[2:17,13]
error_T8 = [(T8_av-T8_lower),-(T8_av-T8_upper)]

T7_av = data[2:17,15]
T7_upper = data[2:17,16]
T7_lower = data[2:17,17]
error_T7 = [(T7_av-T7_lower),-(T7_av-T7_upper)]

T10_av = data[2:17,19]
T10_upper = data[2:17,20]
T10_lower = data[2:17,21]
error_T10 = [(T10_av-T10_lower),-(T10_av-T10_upper)]

T11_av = data[2:17,23]
T11_upper = data[2:17,24]
T11_lower = data[2:17,25]
error_T11 = [(T11_av-T11_lower),-(T11_av-T11_upper)]

#### Growth rate data import ####

grow=genfromtxt('growth_rate_contaminant.csv',delimiter=',')
 
T3 = grow[1:16,1]
T4 = grow[1:16,2]
T8 = grow[1:16,3]
T7 = grow[1:16,4]
T10 = grow[1:16,5]
T11 = grow[1:16,6]


#convert cells/mL to cells/L
T3_av=T3_av*1000
T8_av=T8_av*1000
T10_av=T10_av*1000
T4_av=T4_av*1000
T7_av=T7_av*1000
T11_av=T11_av*1000


figure(10, figsize=(8,6.5))
scatter(day,T3_av, color='#88CCEE', s=60, marker='P', label='RM-3')
errorbar(day,T3_av,yerr=error_T3,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T8_av, color='#228833', s=60, marker='P', label='RM-8')
errorbar(day,T8_av,yerr=error_T8,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T10_av, color='#4477AA', s=60, marker='P', label='RM-10')
errorbar(day,T10_av,yerr=error_T10,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T4_av, color='#DDCC77', s=60, marker='v', label='HM-4')
errorbar(day,T4_av,yerr=error_T4,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T7_av, color='#CC6677', s=60, marker='v', label=' HM-7')
errorbar(day,T7_av,yerr=error_T7,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
scatter(day,T11_av, color='#882255', s=60, marker='v', label= 'HM-11')
errorbar(day,T11_av,yerr=error_T11,xerr=None,fmt='none',color='black',elinewidth=1,capsize=2)
xlabel('Day of Experiment', fontsize=18)
ticklabel_format(axis="y", style="sci", scilimits=(0,0))
xticks(fontsize=16)
yticks(fontsize=16)
ylim(top=5e8)


######### Read in mesocosm data: nutrient information from filters)  #############
##EUKARYOTIC TANKS###

data_nut= genfromtxt('nutrient_cells.csv', delimiter=',')

cell_N = data_nut[1:28,7]
cell_C = data_nut[1:28,9]
cell_P = data_nut[1:28,11]

tank_1_N = array([cell_N[0],cell_N[3],cell_N[9],cell_N[15],cell_N[21]])
tank_2_N = array([cell_N[4],cell_N[10],cell_N[16],cell_N[22]])
tank_5_N = array([cell_N[1],cell_N[5],cell_N[11],cell_N[17],cell_N[23]])
tank_6_N = array([cell_N[6],cell_N[12],cell_N[18],cell_N[24]])
tank_9_N = array([cell_N[7],cell_N[13],cell_N[19],cell_N[25]])
tank_12_N = array([cell_N[2],cell_N[8],cell_N[14],cell_N[20],cell_N[26]])

tank_1_C = array([cell_C[0],cell_C[3],cell_C[9],cell_C[15],cell_C[21]])
tank_2_C = array([cell_C[4],cell_C[10],cell_C[16],cell_C[22]])
tank_5_C = array([cell_C[1],cell_C[5],cell_C[11],cell_C[17],cell_C[23]])
tank_6_C = array([cell_C[6],cell_C[12],cell_C[18],cell_C[24]])
tank_9_C = array([cell_C[7],cell_C[13],cell_C[19],cell_C[25]])
tank_12_C = array([cell_C[2],cell_C[8],cell_C[14],cell_C[20],cell_C[26]])


tank_1_P = array([cell_P[0],cell_P[3],cell_P[9],cell_P[15],cell_P[21]])
tank_2_P = array([cell_P[4],cell_P[10],cell_P[16],cell_P[22]])
tank_5_P = array([cell_P[1],cell_P[5],cell_P[11],cell_P[17],cell_P[23]])
tank_6_P = array([cell_P[6],cell_P[12],cell_P[18],cell_P[24]])
tank_9_P = array([cell_P[7],cell_P[13],cell_P[19],cell_P[25]])
tank_12_P = array([cell_P[2],cell_P[8],cell_P[14],cell_P[20],cell_P[26]])

tank_1_NC = tank_1_N/tank_1_C
tank_2_NC = tank_2_N/tank_2_C
tank_5_NC = tank_5_N/tank_5_C
tank_6_NC = tank_6_N/tank_6_C
tank_9_NC = tank_9_N/tank_9_C
tank_12_NC = tank_12_N/tank_12_C

tank_1_PC = tank_1_P/tank_1_C
tank_2_PC = tank_2_P/tank_2_C
tank_5_PC = tank_5_P/tank_5_C
tank_6_PC = tank_6_P/tank_6_C
tank_9_PC = tank_9_P/tank_9_C
tank_12_PC = tank_12_P/tank_12_C

tank_1_NP = tank_1_N/tank_1_P
tank_2_NP = tank_2_N/tank_2_P
tank_5_NP = tank_5_N/tank_5_P
tank_6_NP = tank_6_N/tank_6_P
tank_9_NP = tank_9_N/tank_9_P
tank_12_NP = tank_12_N/tank_12_P

day= [0,4,7,11,14]
day_1 = [4,7,11,14]

####stoichiometry and day of experiment ######
# 
figure(11, figsize=(8,6.5))
scatter(day,tank_1_NC, color='#88CCEE', s=60,  label='Tank 1')
scatter(day,tank_5_NC, color='#228833', s=60,  label='Tank 5')
scatter(day_1,tank_9_NC, color='#4477AA', s=60,  label='Tank 9')
scatter(day_1,tank_2_NC, color='#DDCC77', s=60,marker='s', label='Tank 2')
scatter(day_1,tank_6_NC, color='#CC6677', s=60,marker='s', label=' Tank 6')
scatter(day,tank_12_NC, color='#882255', s=60,marker='s', label= 'Tank 12')
xticks(fontsize=12)
ylabel('N:C (mol N mol C$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=0.2)

 
figure(12, figsize=(8,6.5))
scatter(day,tank_1_PC, color='#88CCEE', s=60, label='Tank 1')
scatter(day,tank_5_PC, color='#228833', s=60, label='Tank 5')
scatter(day_1,tank_9_PC, color='#4477AA', s=60, label='Tank 9')
scatter(day_1,tank_2_PC, color='#DDCC77', s=60,marker='s', label='Tank 2')
scatter(day_1,tank_6_PC, color='#CC6677', s=60, marker='s',label=' Tank 6')
scatter(day,tank_12_PC, color='#882255', s=60,marker='s', label= 'Tank 12')
title("Eukaryotic Tanks", fontsize=20,x=0.5, y=1.09)
xticks(fontsize=12)
ylabel('P:C (mol P mol C$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=0.02)



figure(13, figsize=(8,6.5))
scatter(day,tank_1_NP, color='#88CCEE', s=60,label='Tank 1')
scatter(day,tank_5_NP, color='#228833', s=60, label='Tank 5')
scatter(day_1,tank_9_NP, color='#4477AA', s=60, label='Tank 9')
scatter(day_1,tank_2_NP, color='#DDCC77', s=60,marker='s', label='Tank 2')
scatter(day_1,tank_6_NP, color='#CC6677', s=60, marker='s',label=' Tank 6')
scatter(day,tank_12_NP, color='#882255', s=60,marker='s', label= 'Tank 12')
xticks(fontsize=12)
ylabel('N:P (mol N mol P$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=70)


######### Mixed population tanks #############
data_nut= genfromtxt('nutrient_cells_contaminant.csv', delimiter=',')

cell_N = data_nut[1:28,5]
cell_C = data_nut[1:28,7]
cell_P = data_nut[1:28,9]

tank_3_N = array([cell_N[0],cell_N[3],cell_N[9],cell_N[15],cell_N[21]])
tank_4_N = array([cell_N[4],cell_N[10],cell_N[16],cell_N[22]])
tank_7_N = array([cell_N[1],cell_N[5],cell_N[11],cell_N[17],cell_N[23]])
tank_8_N = array([cell_N[6],cell_N[12],cell_N[18],cell_N[24]])
tank_10_N = array([cell_N[7],cell_N[13],cell_N[19],cell_N[25]])
tank_11_N = array([cell_N[2],cell_N[8],cell_N[14],cell_N[20],cell_N[26]])


tank_3_C = array([cell_C[0],cell_C[3],cell_C[9],cell_C[15],cell_C[21]])
tank_4_C = array([cell_C[4],cell_C[10],cell_C[16],cell_C[22]])
tank_7_C = array([cell_C[1],cell_C[5],cell_C[11],cell_C[17],cell_C[23]])
tank_8_C = array([cell_C[6],cell_C[12],cell_C[18],cell_C[24]])
tank_10_C = array([cell_C[7],cell_C[13],cell_C[19],cell_C[25]])
tank_11_C = array([cell_C[2],cell_C[8],cell_C[14],cell_C[20],cell_C[26]])


tank_3_P = array([cell_P[0],cell_P[3],cell_P[9],cell_P[15],cell_P[21]])
tank_4_P = array([cell_P[4],cell_P[10],cell_P[16],cell_P[22]])
tank_7_P = array([cell_P[1],cell_P[5],cell_P[11],cell_P[17],cell_P[23]])
tank_8_P = array([cell_P[6],cell_P[12],cell_P[18],cell_P[24]])
tank_10_P = array([cell_P[7],cell_P[13],cell_P[19],cell_P[25]])
tank_11_P = array([cell_P[2],cell_P[8],cell_P[14],cell_P[20],cell_P[26]])

tank_3_NC = tank_3_N/tank_3_C
tank_4_NC = tank_4_N/tank_4_C
tank_7_NC = tank_7_N/tank_7_C
tank_8_NC = tank_8_N/tank_8_C
tank_10_NC = tank_10_N/tank_10_C
tank_11_NC = tank_11_N/tank_11_C

tank_3_PC = tank_3_P/tank_3_C
tank_4_PC = tank_4_P/tank_4_C
tank_7_PC = tank_7_P/tank_7_C
tank_8_PC = tank_8_P/tank_8_C
tank_10_PC = tank_10_P/tank_10_C
tank_11_PC = tank_11_P/tank_11_C

tank_3_NP = tank_3_N/tank_3_P
tank_4_NP = tank_4_N/tank_4_P
tank_7_NP = tank_7_N/tank_7_P
tank_8_NP = tank_8_N/tank_8_P
tank_10_NP = tank_10_N/tank_10_P
tank_11_NP = tank_11_N/tank_11_P

day= [0,4,7,11,14]
day_1 = [4,7,11,14]

####stoichiometry and day of experiment ######
# 
figure(14, figsize=(8,6.5))
scatter(day,tank_3_NC, color='#88CCEE', s=60,marker='P', label='Tank 3')
scatter(day_1,tank_8_NC, color='#228833', s=60, marker='P', label='Tank 8')
scatter(day_1,tank_10_NC, color='#4477AA', s=60, marker='P', label='Tank 10')
scatter(day_1,tank_4_NC, color='#DDCC77', s=60,marker='v', label='Tank 4')
scatter(day,tank_7_NC, color='#CC6677', s=60, marker='v',label=' Tank 7')
scatter(day,tank_11_NC, color='#882255', s=60,marker='v', label= 'Tank 11')
xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
ylabel('N:C (mol N mol C$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=0.2)

figure(15, figsize=(8,6.5))
scatter(day,tank_3_PC, color='#88CCEE', s=60,marker='P', label='Tank 3')
scatter(day_1,tank_8_PC, color='#228833', s=60,marker='P', label='Tank 8')
scatter(day_1,tank_10_PC, color='#4477AA', s=60,marker='P', label='Tank 10')
scatter(day_1,tank_4_PC, color='#DDCC77', s=60,marker='v', label='Tank 4')
scatter(day,tank_7_PC, color='#CC6677', s=60,marker='v', label=' Tank 7')
scatter(day,tank_11_PC, color='#882255', s=60,marker='v', label= 'Tank 11')
xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
title("Mixed Population Tanks", fontsize=20,x=0.5, y=1.09)
ylabel('P:C (mol P mol C$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=0.02)

figure(16, figsize=(8,6.5))
scatter(day,tank_3_NP, color='#88CCEE', s=60,marker='P', label='Tank 3')
scatter(day_1,tank_8_NP, color='#228833', s=60, marker='P',label='Tank 8')
scatter(day_1,tank_10_NP, color='#4477AA', s=60,marker='P', label='Tank 10')
scatter(day_1,tank_4_NP, color='#DDCC77', s=60,marker='v', label='Tank 4')
scatter(day,tank_7_NP, color='#CC6677', s=60,marker='v', label=' Tank 7')
scatter(day,tank_11_NP, color='#882255', s=60,marker='v', label= 'Tank 11')
xlabel('Day of Experiment', fontsize=16)
xticks(fontsize=12)
ylabel('N:P (mol N mol P$^{-1}$)', fontsize=16)
yticks(fontsize=12)
ylim(bottom=0, top=70)




show()
