'''
Created on May 10, 2024

@author: garmin
'''

from pylab import *
from af001_energy_calculation import *
from matplotlib.pyplot import plot, show, legend, figure, scatter, xticks,yticks,ylim,xlim,xlabel,ylabel,title,stackplot,bar,subplots
from matplotlib import *
import matplotlib.patches as mpat
import matplotlib.markers as mar
from sklearn.linear_model import LinearRegression
from numpy import *
from Climitation_max_growth_rate import *
#from sf2 import *

global What_is_limiting

What_is_limiting=0                                          #0: P-limiting  1:N-limiting

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#Function beging here
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
def kkI(I,Pin,Nin,T,G):  
   
    I=I                                                    #change for different irradiances

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if What_is_limiting==0:
    #For P-limiting case
        Pin=Pin#0.002                                           #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=Nin#0.2                                             #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=1.00*10**(-12)/12                                #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    elif What_is_limiting==1:
    #For N-limiting case
        Pin=Pin#0.02                                            #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=Nin#0.05                                            #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=10**(-12)/12                                     #(molC/cell) biomass C per cell (196-18)(from Healey 1985)

    E3=evalue()
    E=E3.E
    
    mu_max=Clim_mu(I)

    C=6
    if G<0:
        Dd=0
    elif G>mu_max:
        Dd=mu_max
    else:
        Dd=G#0.25                                                 #change for different growth rate
    D=Dd/(3600*24)
    Mchl=893.49                                             #(g / mol chlorophyll) mollar mass of chlorophyll
    
    #==============================
    #parameter sets
    #==============================
    m=3.791E-19 #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
    Pmax=0.003205
    OT=0.008634
    Ynphoto_chl=3.561           #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth=4.347E-10           #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=0.5*4.453E-15    #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=2.91679384515998E-15           #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=6213          #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.02816        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=1.5*5.445E-17     #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=25.26/(3.097e16)     
    Cessential=1.518E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
    Ea=96480*0.32#70000  segura, liu2018                              #(J/mol) activation energy for phytoplankton (Li, 1980) cited in (Geider, 1997)
    R=8.3                                   #(J/mol*k) universal gas constant 
    Tref=293                                                #Reference temperature (K) from Healey experiment 20 celsius
    #------------------------------
    #Photosynthesis
    #------------------------------
    Tref=293 
    K = 273                                                 #Reference temperature (K)

    Tt=T+K #arange(Tmin,Tmax+Ttstep,Ttstep)
    U = arange(size(Tt))
    A=Ea/R
    Arr=exp(-A*((1/Tt)-(1/Tref)))                           #arrehenius equation (Geider, 1997) function of temperature 
    
    Parr=ones(size(Arr))                                    #no temperature dependence on photosynthesis
   

    Pchl=Parr*Pmax*(1-exp(-OT*I))                            #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
                                           
    Cnbiosynth = Cnbiosynth/Arr
    Cnrna_variable = Cnrna_variable/Arr
    
    ls=D*Qc                                                 #(molC s-1) Biomass synthesis rate (193-25)
    #================================
    Daverage=0.6*86400
    Cnproteinsynth=24000                                    #(s) Constant for protein synthesis nitrogen calculation (193-25)
    Cndna_variable=1.60e-15                                 #(molN cell-1 s) Constant for variable part of DNA (193-26) 
    Nconst_dnarna=2.78e-15/5/2.8                            #(molN cell-1) Constant part of nitrogen in DNA and RNA (193-26)
    Nconst_other=2.78e-16                                   #(molN cell-1) Other kinds of nitrogen assuming constant (193-33)
    #================================
    
    Molar_mass_DNA_AT_average=307.47                        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97                        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47                        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97                        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506                                          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli                                     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli      #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli      #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=20/7.6                               #(ug/ug) Bremer and Dennis 1996
    RNA_DNA_mass_ratio=17.844/6.5239                        #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli                #(mol mol-1)
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                                               #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    CG=0.563                                                #GC%    [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)                         #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    YdnaC_N=3.5*(1-CG)+2.5*CG                               #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG                              #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

    DNAmb=2.1269                                            #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    Avogadro=6.022*10**23                                   #(molecules mol-1) Avogadro constant
    Pdna_const=DNAmb*2*10**6/Avogadro                       #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio               #(molP cell-1) Constant part of RNA in phosphorus
    #* Make sure to multiply by 2 as they are base PAIRs"
    Ndna_const=Pdna_const/YnucacidP_N                       #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio               #(molN cell-1) Constatn part of RNA in phosphorus
    
    Ynrnaconst_dnarnaconst=1/2                              #(dimensionless) N molar ratio of RNA (constant part) to DNA + RNA (constant part) (193-33) reffering to around p.112 of Biology of Prokyariotes
    Yndnaconst_dnarnaconst=1-Ynrnaconst_dnarnaconst         #(dimensionless) N molar ratio of DNA (constant part) to DNA + RNA (constant part)  (193-33) refering to around p.112 of Biology of Prokyariotes
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Chl=((1+E)*ls+m)/Pchl                                   #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl                                  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nproteinsynth=(Nphoto+Nconst_protein)*D*Cnproteinsynth/(1-D*Cnproteinsynth)                       #(molN cell-1) protein synthesis related protein in N (193-25)
    Nbiosynth=D*Cnbiosynth                                  #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth                #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable                 #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=Ndna_const*Dd/1.2*(18.3-7.6)/7.6          #(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996
    Ndna_variable=0*Dd                                      #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
    Nchl=Chl*4/55                                           #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl                          #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N                 #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N                 #(molP cell-1) variable part of phosphorus in DNA (193-26)
    
    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max               #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)
                                                                                                    #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl                          #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const                #(molP cell-1) total phosphorus in the cell without storage
    
    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(Tt))
    X=zeros(size(Tt))
    Qn_test=zeros(size(Tt))
    Qp_test=copy(X)
    Qp=copy(X)
    Qn=copy(X)
    Pstore=copy(X)
    #=================================
    #Population calculation
    #=================================
    Xn_max=Nin/Qn_min
    Xp_max=Pin/Qp_min
    for i in U:
        if Xn_max[i]>Xp_max[i]:
            X[i]=Xp_max[i]
            Qp[i]=Qp_min[i]
            Qn_test[i]=Nin/X[i]
            if Qn_test[i]<Qn_max[i]:
                Qn[i]=Qn_test[i]
                Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]    #(molN cell-1) Nitrogen storage in the cell
            else:
                Qn[i]=Qn_max[i]
                Nstore[i]=Nstore_max
        else:
            X[i]=Xn_max[i]
            Qn[i]=Qn_min[i]
            Qp_test[i]=Pin/X[i]
            if Qp_test[i]<Qp_max:
                Qp[i]=Qp_test[i]
            else:
                Qp[i]=Qp_max
            Pstore[i]=Qp[i]-Pconst_other-Pthylakoid[i]-Prna_variable[i]-Prna_const-Pdna_variable-Pdna_const         #(molP cell-1) Stored phosphorus in the cell
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc                            #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc*14*10**6/(12*10**3)          #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc*30.97*10**6/(12*10**3)       #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp*14*10**6/(30.97*10**6)       #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc                              #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49                                 #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6    #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 3 (unit adjustment)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nunit=1/Qc                                  #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
    Punit=1/Qc                                  #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
    Numbertoarray=ones(size(Tt))                #(dimensionless) Number to array converter
    
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    CNprotein=4.49                              #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    
    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein                     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein               #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein     #(molC cell-1) carbon in other protein assumed constant (195-16)
 
    #----------------------
    #C chlorophyll
    #----------------------
    Cchl=Chl                                    #(molC cell-1) carbon in chlorophyll (195-16)
    
    #----------------------
    #C DNA RNA
    #----------------------
    Crna_const=Nrna_const*YrnaC_N               #(molC cell-1) carbon in variable part of RNA (195-16)
    Crna_variable=Nrna_variable*YrnaC_N         #(molC cell-1) carbon in variable part of RNA (195-16)
    
    Cdna_const=Ndna_const*YdnaC_N               #(molC cell-1) carbon in constant part of DNA (195-16)
    Cdna_variable=Ndna_variable*YdnaC_N         #(molC cell-1) carbon in variable part of DNA (195-16)
    

    Cnstore=Nstore*YcyanoC_N                    #(molC cell-1) carbon in nitrogen storage (cyanophycin)
    CthylakoidPG=Pthylakoid*YpgC_P              #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    
    #---------------------------------------------------
    #C other: Here revised to include Nstore reduction
    #---------------------------------------------------

    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
    Cother=Cother_with_full_Nstore            

    Cother[Cother<0]=0
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc                            #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc#*14*10**6/(12*10**3)          #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc#*30.97*10**6/(12*10**3)       #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp                              #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc                              #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49                                 #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6    #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #=======================================
    #C for plot in
    #=======================================
    percentorratio=100                          #100: percent, 1:ratio
    Cphoto_plot=Cphoto/Qc*percentorratio           
    Cbiosynth_plot=Cbiosynth/Qc*percentorratio
    Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numbertoarray
    Cchl_plot=Cchl/Qc*percentorratio
    Crna_const_plot=Crna_const/Qc*percentorratio*Numbertoarray
    Crna_variable_plot=Crna_variable/Qc*percentorratio
    Cdna_const_plot=Cdna_const/Qc*percentorratio*Numbertoarray
    Cdna_variable_plot=Cdna_variable/Qc*percentorratio

    Cother_plot=Cother/Qc*percentorratio
    Cessential_plot=Cessential/Qc*percentorratio*Numbertoarray
    Cnstore_plot=Cnstore/Qc*percentorratio
    CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio
    
    if Cphoto_plot>100:
        Cphoto_plot=100-Cbiosynth_plot-Cconst_protein_plot-Cchl_plot-Crna_const_plot-Crna_variable_plot-Cdna_const_plot-Cdna_variable_plot-Cother_plot-Cessential_plot-Cnstore_plot-CthylakoidPG_plot
    
    Nphoto_plot=Nphoto*Nunit    #(ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
    Nbiosynth_plot=Nbiosynth*Nunit      #(ug N/ mgC) biosynthesis related protein in N (193-37)
    Nconst_protein_plot=Nconst_protein*Nunit*Numbertoarray    #(ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
    Nchl_plot=Nchl*Nunit        #(ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)

    Nrna_variable_plot=Nrna_variable*Nunit      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    Ndna_variable_plot=Ndna_variable*Nunit      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    
    Ndna_const_plot=Ndna_const*Nunit*Numbertoarray    #(ug N/ mgC) Nitrogen in constant part of DNA
    Nrna_const_plot=Nrna_const*Nunit*Numbertoarray    #(ug N/ mgC) Nitrogen in constant part of RNA
    
    Nstore_plot=Nstore*Nunit      #(ug N/ mgC) Nitrogen in storage
    #=======================================
    #P for plot ***(Mainly from 193-33)***
    #=======================================
    Prna_variable_plot=Prna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of RNA (193-37)
    Pdna_variable_plot=Pdna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of DNA (193-37)
    
    Pthylakoid_plot=Pthylakoid*Punit     #(ug P/ mgC) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)(193-33)
    Pconst_other_plot=Pconst_other*Punit*Numbertoarray       #(ug P/ mgC) PHosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)
    Pdna_const_plot=Pdna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of DNA
    Prna_const_plot=Prna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of RNA
    
    Pstore_plot=Pstore*Punit
    

    Nbio=Nbiosynth_plot+Nrna_variable_plot+Nrna_const_plot
    Nphoto=Nphoto_plot+Nchl_plot
    Ness=Nconst_protein_plot+Ndna_const_plot+Ndna_variable_plot
    Nsto=Nstore_plot
    
    Pbio=Prna_variable_plot+Prna_const_plot
    Pphoto=Pthylakoid_plot
    Pess=Pconst_other_plot+Pdna_const_plot+Pdna_variable_plot
    Psto=Pstore_plot


    return Nbio,Nphoto,Ness,Nsto,Pbio,Pphoto,Pess,Psto

env= genfromtxt('model_input.csv',delimiter=',')
time=env[1:2161,1]
 
I1=(1/(10.8*5))*env[1:2161,3]
I5=(1/(10.8*5))*env[1:2161,8]
I9=(1/(10.8*5))*env[1:2161,13]
I2=(1/(10.8*5))*env[1:2161,18]
I6=(1/(10.8*5))*env[1:2161,23]
I12=(1/(10.8*5))*env[1:2161,28]

#convert ug/l to mol/m3
Pin_1=env[1:2161,5]*(1000/(30.97*1e6))
Pin_5=env[1:2161,10]*(1000/(30.97*1e6))
Pin_9=env[1:2161,15]*(1000/(30.97*1e6))
Pin_2=env[1:2161,20]*(1000/(30.97*1e6))
Pin_6=env[1:2161,25]*(1000/(30.97*1e6))
Pin_12=env[1:2161,30]*(1000/(30.97*1e6))


Nin_1=env[1:2161,4]*(1000/(14*1e6))
Nin_5=env[1:2161,9]*(1000/(14*1e6))
Nin_9=env[1:2161,14]*(1000/(14*1e6))
Nin_2=env[1:2161,19]*(1000/(14*1e6))
Nin_6=env[1:2161,24]*(1000/(14*1e6))
Nin_12=env[1:2161,29]*(1000/(14*1e6))

G1 = env[1:2161,6]
G5 = env[1:2161,11]
G9 = env[1:2161,16]
G2 = env[1:2161,21]
G6 = env[1:2161,26]
G12 = env[1:2161,31]

T1=env[1:2161,2]
T5=env[1:2161,7]
T9=env[1:2161,12]
T2=env[1:2161,17]
T6=env[1:2161,22]
T12=env[1:2161,27]

def use_fnx(T,I,Pin,Nin,G):
    def day(a):
        return a[0:144],a[144:288],a[288:432],a[432:576],a[576:720],a[720:864],a[864:1008],a[1008:1152],a[1152:1296],a[1296:1440],a[1440:1584],a[1584:1728],a[1728:1872],a[1872:2016],a[2016:]
 
    def morning_samp(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o):
        return a[42:55],b[42:55],c[42:55],d[42:55],e[42:55],f[42:55],g[42:55],h[42:55],i[42:55],j[42:55],k[42:55],l[42:55],m[42:55],n[42:55],o[42:55]
 
    def avg(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o):
        return array([sum(a)/len(a),sum(b)/len(b),sum(c)/len(c),sum(d)/len(d),sum(e)/len(e),sum(f)/len(f),sum(g)/len(g),sum(h)/len(h),sum(i)/len(i),sum(j)/len(j),sum(k)/len(k),sum(l)/len(l),sum(m)/len(m),sum(n)/len(n),sum(o)/len(o)])
 
    d0_T1,d1_T1,d2_T1,d3_T1,d4_T1,d5_T1,d6_T1,d7_T1,d8_T1,d9_T1,d10_T1,d11_T1,d12_T1,d13_T1,d14_T1=day(T)
    d0_I1,d1_I1,d2_I1,d3_I1,d4_I1,d5_I1,d6_I1,d7_I1,d8_I1,d9_I1,d10_I1,d11_I1,d12_I1,d13_I1,d14_I1=day(I)
    d0_Pin_1,d1_Pin_1,d2_Pin_1,d3_Pin_1,d4_Pin_1,d5_Pin_1,d6_Pin_1,d7_Pin_1,d8_Pin_1,d9_Pin_1,d10_Pin_1,d11_Pin_1,d12_Pin_1,d13_Pin_1,d14_Pin_1=day(Pin)
    d0_Nin_1,d1_Nin_1,d2_Nin_1,d3_Nin_1,d4_Nin_1,d5_Nin_1,d6_Nin_1,d7_Nin_1,d8_Nin_1,d9_Nin_1,d10_Nin_1,d11_Nin_1,d12_Nin_1,d13_Nin_1,d14_Nin_1=day(Nin)
    d0_G1,d1_G1,d2_G1,d3_G1,d4_G1,d5_G1,d6_G1,d7_G1,d8_G1,d9_G1,d10_G1,d11_G1,d12_G1,d13_G1,d14_G1=day(G)
    
    m0_T1,m1_T1,m2_T1,m3_T1,m4_T1,m5_T1,m6_T1,m7_T1,m8_T1,m9_T1,m10_T1,m11_T1,m12_T1,m13_T1,m14_T1=morning_samp(d0_T1,d1_T1,d2_T1,d3_T1,d4_T1,d5_T1,d6_T1,d7_T1,d8_T1,d9_T1,d10_T1,d11_T1,d12_T1,d13_T1,d14_T1)
    m0_I1,m1_I1,m2_I1,m3_I1,m4_I1,m5_I1,m6_I1,m7_I1,m8_I1,m9_I1,m10_I1,m11_I1,m12_I1,m13_I1,m14_I1=morning_samp(d0_I1,d1_I1,d2_I1,d3_I1,d4_I1,d5_I1,d6_I1,d7_I1,d8_I1,d9_I1,d10_I1,d11_I1,d12_I1,d13_I1,d14_I1)
    m0_Pin_1,m1_Pin_1,m2_Pin_1,m3_Pin_1,m4_Pin_1,m5_Pin_1,m6_Pin_1,m7_Pin_1,m8_Pin_1,m9_Pin_1,m10_Pin_1,m11_Pin_1,m12_Pin_1,m13_Pin_1,m14_Pin_1=morning_samp(d0_Pin_1,d1_Pin_1,d2_Pin_1,d3_Pin_1,d4_Pin_1,d5_Pin_1,d6_Pin_1,d7_Pin_1,d8_Pin_1,d9_Pin_1,d10_Pin_1,d11_Pin_1,d12_Pin_1,d13_Pin_1,d14_Pin_1)
    m0_Nin_1,m1_Nin_1,m2_Nin_1,m3_Nin_1,m4_Nin_1,m5_Nin_1,m6_Nin_1,m7_Nin_1,m8_Nin_1,m9_Nin_1,m10_Nin_1,m11_Nin_1,m12_Nin_1,m13_Nin_1,m14_Nin_1=morning_samp(d0_Nin_1,d1_Nin_1,d2_Nin_1,d3_Nin_1,d4_Nin_1,d5_Nin_1,d6_Nin_1,d7_Nin_1,d8_Nin_1,d9_Nin_1,d10_Nin_1,d11_Nin_1,d12_Nin_1,d13_Nin_1,d14_Nin_1)
    m0_G1,m1_G1,m2_G1,m3_G1,m4_G1,m5_G1,m6_G1,m7_G1,m8_G1,m9_G1,m10_G1,m11_G1,m12_G1,m13_G1,m14_G1=morning_samp(d0_G1,d1_G1,d2_G1,d3_G1,d4_G1,d5_G1,d6_G1,d7_G1,d8_G1,d9_G1,d10_G1,d11_G1,d12_G1,d13_G1,d14_G1)
    
    avg_T_1=avg(m0_T1,m1_T1,m2_T1,m3_T1,m4_T1,m5_T1,m6_T1,m7_T1,m8_T1,m9_T1,m10_T1,m11_T1,m12_T1,m13_T1,m14_T1)
    avg_I_1=avg(m0_I1,m1_I1,m2_I1,m3_I1,m4_I1,m5_I1,m6_I1,m7_I1,m8_I1,m9_I1,m10_I1,m11_I1,m12_I1,m13_I1,m14_I1)
    avg_Pin_1=avg(m0_Pin_1,m1_Pin_1,m2_Pin_1,m3_Pin_1,m4_Pin_1,m5_Pin_1,m6_Pin_1,m7_Pin_1,m8_Pin_1,m9_Pin_1,m10_Pin_1,m11_Pin_1,m12_Pin_1,m13_Pin_1,m14_Pin_1)
    avg_Nin_1=avg(m0_Nin_1,m1_Nin_1,m2_Nin_1,m3_Nin_1,m4_Nin_1,m5_Nin_1,m6_Nin_1,m7_Nin_1,m8_Nin_1,m9_Nin_1,m10_Nin_1,m11_Nin_1,m12_Nin_1,m13_Nin_1,m14_Nin_1)
    avg_G_1=avg(m0_G1,m1_G1,m2_G1,m3_G1,m4_G1,m5_G1,m6_G1,m7_G1,m8_G1,m9_G1,m10_G1,m11_G1,m12_G1,m13_G1,m14_G1)
    
    return avg_T_1,avg_I_1,avg_Pin_1,avg_Nin_1,avg_G_1


avg_T_1,avg_I_1,avg_Pin_1,avg_Nin_1,avg_G_1=use_fnx(T1,I1,Pin_1,Nin_1,G1)
avg_T_5,avg_I_5,avg_Pin_5,avg_Nin_5,avg_G_5=use_fnx(T5,I5,Pin_5,Nin_5,G5)
avg_T_9,avg_I_9,avg_Pin_9,avg_Nin_9,avg_G_9=use_fnx(T9,I9,Pin_9,Nin_9,G9)
avg_T_2,avg_I_2,avg_Pin_2,avg_Nin_2,avg_G_2=use_fnx(T2,I2,Pin_2,Nin_2,G2)
avg_T_6,avg_I_6,avg_Pin_6,avg_Nin_6,avg_G_6=use_fnx(T6,I6,Pin_6,Nin_6,G6)
avg_T_12,avg_I_12,avg_Pin_12,avg_Nin_12,avg_G_12=use_fnx(T12,I12,Pin_12,Nin_12,G12)



def mod_macro_stoich(avg_I,avg_Pin,avg_Nin,avg_T,avg_G):
    loop=arange(0,15,1)
    Nb=[]
    Np=[]
    Ne=[]
    Ns=[]
    Pb=[]
    Pp=[]
    Pe=[]
    Ps=[]
    
    for i in loop:
        nb,np,ne,ns,pb,pp,pe,ps=kkI(avg_I[i],avg_Pin[i],avg_Nin[i],avg_T[i],avg_G[i])
        #print(cb)
        Nb.append(nb.item())
        Np.append(np.item())
        Ne.append(ne.item())
        Ns.append(ns.item())
        Pb.append(pb.item())
        Pp.append(pp.item())
        Pe.append(pe.item())
        Ps.append(ps.item())
    
    Nb=array(Nb)
    Np=array(Np)
    Ne=array(Ne)
    Ns=array(Ns)
    Pb=array(Pb)
    Pp=array(Pp)
    Pe=array(Pe)
    Ps=array(Ps)

    
    return Nb,Np,Ne,Ns,Pb,Pp,Pe,Ps


Nb1,Np1,Ne1,Ns1,Pb1,Pp1,Pe1,Ps1=mod_macro_stoich(avg_I_1,avg_Pin_1,avg_Nin_1,avg_T_1,avg_G_1)
Nb5,Np5,Ne5,Ns5,Pb5,Pp5,Pe5,Ps5=mod_macro_stoich(avg_I_5,avg_Pin_5,avg_Nin_5,avg_T_5,avg_G_5)
Nb9,Np9,Ne9,Ns9,Pb9,Pp9,Pe9,Ps9=mod_macro_stoich(avg_I_9,avg_Pin_9,avg_Nin_9,avg_T_9,avg_G_9)
Nb2,Np2,Ne2,Ns2,Pb2,Pp2,Pe2,Ps2=mod_macro_stoich(avg_I_2,avg_Pin_2,avg_Nin_2,avg_T_2,avg_G_2)
Nb6,Np6,Ne6,Ns6,Pb6,Pp6,Pe6,Ps6=mod_macro_stoich(avg_I_6,avg_Pin_6,avg_Nin_6,avg_T_6,avg_G_6)
Nb12,Np12,Ne12,Ns12,Pb12,Pp12,Pe12,Ps12=mod_macro_stoich(avg_I_12,avg_Pin_12,avg_Nin_12,avg_T_12,avg_G_12)



######### Read in mesocosm data: nutrient information from filters)  #############
# ######### Read in mesocosm data: nutrient information from filters)  #############
data_nut= genfromtxt('nutrient_cells_contaminant.csv', delimiter=',')
 
cell_N = data_nut[1:28,5]
cell_C = data_nut[1:28,7]
cell_P = data_nut[1:28,9]
 
tank_1_N = array([cell_N[0],cell_N[3],cell_N[9],cell_N[15],cell_N[21]])
tank_2_N = array([cell_N[4],cell_N[10],cell_N[16],cell_N[22]])
tank_6_N = array([cell_N[1],cell_N[5],cell_N[11],cell_N[17],cell_N[23]])
tank_5_N = array([cell_N[6],cell_N[12],cell_N[18],cell_N[24]])
tank_9_N = array([cell_N[7],cell_N[13],cell_N[19],cell_N[25]])
tank_12_N = array([cell_N[2],cell_N[8],cell_N[14],cell_N[20],cell_N[26]])
 
tank_1_C = array([cell_C[0],cell_C[3],cell_C[9],cell_C[15],cell_C[21]])
tank_2_C = array([cell_C[4],cell_C[10],cell_C[16],cell_C[22]])
tank_6_C  = array([cell_C[1],cell_C[5],cell_C[11],cell_C[17],cell_C[23]])
tank_5_C= array([cell_C[6],cell_C[12],cell_C[18],cell_C[24]])
tank_9_C = array([cell_C[7],cell_C[13],cell_C[19],cell_C[25]])
tank_12_C = array([cell_C[2],cell_C[8],cell_C[14],cell_C[20],cell_C[26]])
 
 
tank_1_P = array([cell_P[0],cell_P[3],cell_P[9],cell_P[15],cell_P[21]])
tank_2_P = array([cell_P[4],cell_P[10],cell_P[16],cell_P[22]])
tank_6_P = array([cell_P[1],cell_P[5],cell_P[11],cell_P[17],cell_P[23]])
tank_5_P  = array([cell_P[6],cell_P[12],cell_P[18],cell_P[24]])
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
 


def NC_allo_figure(Nb,Np,Ne,Ns,NC,sam_num,name,tank):
    days=arange(0,15,1)
    ind = arange(15)
    width = 0.75 
    figure(figsize=(8, 6.5))

    
    bar(days, Nb, width=width, label='Bio', color='#44AA99', bottom=Np+Ns+Ne)
    bar(days, Np, width=width, label='Pho', color='#CC6677', bottom=Ns+Ne)
    bar(days, Ns, width=width, label='N.Sto', color='#882255', bottom=Ne)
    bar(days, Ne, width=width, label='Ess', color='#4B0082')
    

    xticks(days,days,fontsize=16)

    xlabel("Day of Experiment",fontsize=18)
    ylabel('N:C (mol N mol C$^{-1}$)',fontsize=18)
    yticks(fontsize=18)
    ylim(0,0.35)
    
    if sam_num==5:
        sday=[0,4,7,11,14]
    else:
        sday=[4,7,11,14]
    
    scatter(sday,NC,color='#999933',edgecolor='black',s=60)

    title(tank,fontsize=20)
    #sf(name,500)
    
  
NC_allo_figure(Nb1,Np1,Ne1,Ns1,tank_1_NC,5,'Tank3_NCallo','RM-3')
NC_allo_figure(Nb5,Np5,Ne5,Ns5,tank_5_NC,4,'Tank8_NCallo','RM-8')
NC_allo_figure(Nb9,Np9,Ne9,Ns9,tank_9_NC,4,'Tank10_NCallo','RM-10')
NC_allo_figure(Nb2,Np2,Ne2,Ns2,tank_2_NC,4,'Tank4_NCallo','HM-4')
NC_allo_figure(Nb6,Np6,Ne6,Ns6,tank_6_NC,5,'Tank7_NCallo','HM-7')
NC_allo_figure(Nb12,Np12,Ne12,Ns12,tank_12_NC,5,'Tank11_NCallo','HM-11')

def PC_allo_figure(Pb,Pp,Pe,Ps,PC,sam_num,name,tank):
    days=arange(0,15,1)
    ind = arange(15)
    width = 0.75 
    figure(figsize=(8, 6.5))

    
    bar(days, Pb, width=width, label='Bio', color='#44AA99', bottom=Pp+Ps+Pe)
    bar(days, Pp, width=width, label='Pho', color='#CC6677', bottom=Ps+Pe)
    bar(days, Ps, width=width, label='P.Sto', color='#88CCEE', bottom=Pe)
    bar(days, Pe, width=width, label='Ess', color='#4B0082')
    

    xticks(days,days,fontsize=16)

    xlabel("Day of Experiment",fontsize=18)
    ylabel('P:C (mol P mol C$^{-1}$)',fontsize=18)
    yt=arange(0,0.025,0.010)
    yticks(ticks=yt,fontsize=18)
    ylim(0,0.02)
    
    if sam_num==5:
        sday=[0,4,7,11,14]
    else:
        sday=[4,7,11,14]
    
    scatter(sday,PC,color='#88CCEE',edgecolor='black',s=60)

    title(tank,fontsize=20)
    #sf(name,500)

PC_allo_figure(Pb1,Pp1,Pe1,Ps1,tank_1_PC,5,'Tank3_PCallo','RM-3')
PC_allo_figure(Pb5,Pp5,Pe5,Ps5,tank_5_PC,4,'Tank8_PCallo','RM-8')
PC_allo_figure(Pb9,Pp9,Pe9,Ps9,tank_9_PC,4,'Tank10_PCallo','RM-10')
PC_allo_figure(Pb2,Pp2,Pe2,Ps2,tank_2_PC,4,'Tank4_PCallo','HM-4')
PC_allo_figure(Pb6,Pp6,Pe6,Ps6,tank_6_PC,5,'Tank7_PCallo','HM-7')
PC_allo_figure(Pb12,Pp12,Pe12,Ps12,tank_12_PC,5,'Tank11_PCallo','HM-11')


show()
