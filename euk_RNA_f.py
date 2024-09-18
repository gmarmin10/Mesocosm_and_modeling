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
    m=3*3.791E-19 #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
    Pmax=0.003205
    OT=0.008634
    Ynphoto_chl=3.561           #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth=0.5*4.347E-10           #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=0.1*4.453E-15    #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=2.91679384515998E-15           #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=2.5*6213          #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.02816        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=5*5.445E-17     #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=5*25.26/(3.097e16)     
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
    
    Ctot=Cphoto_plot+Cbiosynth_plot+Cconst_protein_plot+Cchl_plot+Crna_const_plot+Crna_variable_plot+Cdna_const_plot+Cdna_variable_plot+Cessential_plot+Cnstore_plot+CthylakoidPG_plot
    
    if Ctot>100:
        Cphoto_plot=100-Cbiosynth_plot-Cconst_protein_plot-Cchl_plot-Crna_const_plot-Crna_variable_plot-Cdna_const_plot-Cdna_variable_plot-Cother_plot-Cessential_plot-Cnstore_plot-CthylakoidPG_plot
   

    Cbio=Crna_const_plot+Crna_variable_plot+Cbiosynth_plot
    Cphoto=Cphoto_plot+CthylakoidPG_plot+Cchl_plot
    Cess=Cessential_plot+Cconst_protein_plot+Cdna_const_plot+Cdna_variable_plot
    Csto=Cother_plot
    Cnsto=Cnstore_plot
        
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
    RNA_C=Molar_mass_RNA_AT_average*CG+Molar_mass_RNA_CG_average*(1-CG)
    CRNA=Crna_const_plot+Crna_variable_plot
    CRNA=(Crna_const+Crna_variable)#*(1/10.7)*(RNA_C)
    NRNA=Nrna_variable_plot+Nrna_const_plot
    
    PRNA=Prna_variable_plot+Prna_const_plot

    
    return CRNA,NRNA,PRNA,Cchl

def kkI2(I,Pin,Nin,T,G):  
   
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
    Pconst_other=5.445E-17     #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
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
    
    Ctot=Cphoto_plot+Cbiosynth_plot+Cconst_protein_plot+Cchl_plot+Crna_const_plot+Crna_variable_plot+Cdna_const_plot+Cdna_variable_plot+Cessential_plot+Cnstore_plot+CthylakoidPG_plot
    
    if Ctot>100:
        Cphoto_plot=100-Cbiosynth_plot-Cconst_protein_plot-Cchl_plot-Crna_const_plot-Crna_variable_plot-Cdna_const_plot-Cdna_variable_plot-Cother_plot-Cessential_plot-Cnstore_plot-CthylakoidPG_plot
   

    Cbio=Crna_const_plot+Crna_variable_plot+Cbiosynth_plot
    Cphoto=Cphoto_plot+CthylakoidPG_plot+Cchl_plot
    Cess=Cessential_plot+Cconst_protein_plot+Cdna_const_plot+Cdna_variable_plot
    Csto=Cother_plot
    Cnsto=Cnstore_plot
        
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
    RNA_C=Molar_mass_RNA_AT_average*CG+Molar_mass_RNA_CG_average*(1-CG)
    CRNA=Crna_const_plot+Crna_variable_plot
    CRNA=(Crna_const+Crna_variable)#*(1/10.7)*(RNA_C)
    NRNA=Nrna_variable_plot+Nrna_const_plot
    
    PRNA=Prna_variable_plot+Prna_const_plot

    
    return CRNA,NRNA,PRNA,Cchl

env= genfromtxt('model_input_02.csv',delimiter=',')
time=env[1:2161,1]
 
I1=(1/(10.8*5))*env[1:2161,3]
I5=(1/(10.8*5))*env[1:2161,8]
I9=(1/(10.8*5))*env[1:2161,13]
I2=(1/(10.8*5))*env[1:2161,18]
I6=(1/(10.8*5))*env[1:2161,23]
I12=(1/(10.8*5))*env[1:2161,28]

#convert ug/l to mol/m3
Pin_1=env[1:2161,5]*(1000/(1e6*30.97))
Pin_5=env[1:2161,10]*(1000/(1e6*30.97))
Pin_9=env[1:2161,15]*(1000/(1e6*30.97))
Pin_2=env[1:2161,20]*(1000/(1e6*30.97))
Pin_6=env[1:2161,25]*(1000/(1e6*30.97))
Pin_12=env[1:2161,30]*(1000/(1e6*30.97))


Nin_1=env[1:2161,4]*(1000/(1e6*14))
Nin_5=env[1:2161,9]*(1000/(1e6*14))
Nin_9=env[1:2161,14]*(1000/(1e6*14))
Nin_2=env[1:2161,19]*(1000/(1e6*14))
Nin_6=env[1:2161,24]*(1000/(1e6*14))
Nin_12=env[1:2161,29]*(1000/(1e6*14))

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
    CRNA=[]
    NRNA=[]
    PRNA=[]
    Cchl=[]

    for i in loop:
        crna,nrna,prna,cchl=kkI(avg_I[i],avg_Pin[i],avg_Nin[i],avg_T[i],avg_G[i])
        #print(cb)
        CRNA.append(crna.item())
        NRNA.append(nrna.item())
        PRNA.append(prna.item())
        Cchl.append(cchl.item())

    
    CRNA=array(CRNA)
    NRNA=array(NRNA)
    PRNA=array(PRNA)
    Cchl=array(Cchl)

    
    return CRNA,NRNA,PRNA,Cchl

CRNA1,NRNA1,PRNA1,Cchl1=mod_macro_stoich(avg_I_1,avg_Pin_1,avg_Nin_1,avg_T_1,avg_G_1)
CRNA5,NRNA5,PRNA5,Cchl5=mod_macro_stoich(avg_I_5,avg_Pin_5,avg_Nin_5,avg_T_5,avg_G_5)
CRNA9,NRNA9,PRNA9,Cchl9=mod_macro_stoich(avg_I_9,avg_Pin_9,avg_Nin_9,avg_T_9,avg_G_9)
CRNA2,NRNA2,PRNA2,Cchl2=mod_macro_stoich(avg_I_2,avg_Pin_2,avg_Nin_2,avg_T_2,avg_G_2)
CRNA6,NRNA6,PRNA6,Cchl6=mod_macro_stoich(avg_I_6,avg_Pin_6,avg_Nin_6,avg_T_6,avg_G_6)
CRNA12,NRNA12,PRNA12,Cchl12=mod_macro_stoich(avg_I_12,avg_Pin_12,avg_Nin_12,avg_T_12,avg_G_12)

def mod_macro_stoich2(avg_I,avg_Pin,avg_Nin,avg_T,avg_G):
    loop=arange(0,15,1)
    CRNA=[]
    NRNA=[]
    PRNA=[]
    Cchl=[]

    for i in loop:
        crna,nrna,prna,cchl=kkI2(avg_I[i],avg_Pin[i],avg_Nin[i],avg_T[i],avg_G[i])
        #print(cb)
        CRNA.append(crna.item())
        NRNA.append(nrna.item())
        PRNA.append(prna.item())
        Cchl.append(cchl.item())

    
    CRNA=array(CRNA)
    NRNA=array(NRNA)
    PRNA=array(PRNA)
    Cchl=array(Cchl)

    
    return CRNA,NRNA,PRNA,Cchl

env= genfromtxt('model_input.csv',delimiter=',')
time=env[1:2161,1]

I3=(1/(10.8*5))*env[1:2161,3]
I8=(1/(10.8*5))*env[1:2161,8]
I10=(1/(10.8*5))*env[1:2161,13]
I4=(1/(10.8*5))*env[1:2161,18]
I7=(1/(10.8*5))*env[1:2161,23]
I11=(1/(10.8*5))*env[1:2161,28]

#convert ug/l to mol/m3
Pin_3=env[1:2161,5]*(1000/(1e6*30.97))
Pin_8=env[1:2161,10]*(1000/(1e6*30.97))
Pin_10=env[1:2161,15]*(1000/(1e6*30.97))
Pin_4=env[1:2161,20]*(1000/(1e6*30.97))
Pin_7=env[1:2161,25]*(1000/(1e6*30.97))
Pin_11=env[1:2161,30]*(1000/(1e6*30.97))


Nin_3=env[1:2161,4]*(1000/(1e6*14))
Nin_8=env[1:2161,9]*(1000/(1e6*14))
Nin_10=env[1:2161,14]*(1000/(1e6*14))
Nin_4=env[1:2161,19]*(1000/(1e6*14))
Nin_7=env[1:2161,24]*(1000/(1e6*14))
Nin_11=env[1:2161,29]*(1000/(1e6*14))

G3 = env[1:2161,6]
G8 = env[1:2161,11]
G10 = env[1:2161,16]
G4= env[1:2161,21]
G7 = env[1:2161,26]
G11 = env[1:2161,31]

T3=env[1:2161,2]
T8=env[1:2161,7]
T10=env[1:2161,12]
T4=env[1:2161,17]
T7=env[1:2161,22]
T11=env[1:2161,27]




avg_T_3,avg_I_3,avg_Pin_3,avg_Nin_3,avg_G_3=use_fnx(T3,I3,Pin_3,Nin_3,G3)
avg_T_8,avg_I_8,avg_Pin_8,avg_Nin_8,avg_G_8=use_fnx(T8,I8,Pin_8,Nin_8,G8)
avg_T_10,avg_I_10,avg_Pin_10,avg_Nin_10,avg_G_10=use_fnx(T10,I10,Pin_10,Nin_10,G10)
avg_T_4,avg_I_4,avg_Pin_4,avg_Nin_4,avg_G_4=use_fnx(T4,I4,Pin_4,Nin_4,G4)
avg_T_7,avg_I_7,avg_Pin_7,avg_Nin_7,avg_G_7=use_fnx(T7,I7,Pin_7,Nin_7,G7)
avg_T_11,avg_I_11,avg_Pin_11,avg_Nin_11,avg_G_11=use_fnx(T11,I11,Pin_11,Nin_11,G11)


CRNA3,NRNA3,PRNA3,Cchl3=mod_macro_stoich2(avg_I_3,avg_Pin_3,avg_Nin_3,avg_T_3,avg_G_3)
CRNA8,NRNA8,PRNA8,Cchl8=mod_macro_stoich2(avg_I_8,avg_Pin_8,avg_Nin_8,avg_T_8,avg_G_8)
CRNA10,NRNA10,PRNA10,Cchl10=mod_macro_stoich2(avg_I_10,avg_Pin_10,avg_Nin_10,avg_T_10,avg_G_10)
CRNA4,NRNA4,PRNA4,Cchl4=mod_macro_stoich2(avg_I_4,avg_Pin_4,avg_Nin_4,avg_T_4,avg_G_4)
CRNA7,NRNA7,PRNA7,Cchl7=mod_macro_stoich2(avg_I_7,avg_Pin_7,avg_Nin_7,avg_T_7,avg_G_7)
CRNA11,NRNA11,PRNA11,Cchl11=mod_macro_stoich2(avg_I_11,avg_Pin_11,avg_Nin_11,avg_T_11,avg_G_11)



data=genfromtxt('RNA.csv',delimiter=',') #ng RNA per microliter

tank1_rna_d=data[1:,1]
tank1_rna=array([tank1_rna_d[0],tank1_rna_d[2],tank1_rna_d[4],tank1_rna_d[6]])

tank5_rna_d=data[1:,2]
tank5_rna=array([tank5_rna_d[0],tank5_rna_d[2],tank5_rna_d[4],tank5_rna_d[6]])

tank9_rna_d=data[1:,3]
tank9_rna=array([tank9_rna_d[0],tank9_rna_d[2],tank9_rna_d[4],tank9_rna_d[6]])

tank2_rna_d=data[1:,4]
tank2_rna=array([tank2_rna_d[0],tank2_rna_d[2],tank2_rna_d[4],tank2_rna_d[6]])

tank6_rna_d=data[1:,5]
tank6_rna=array([tank6_rna_d[0],tank6_rna_d[2],tank6_rna_d[4],tank6_rna_d[6]])

tank12_rna_d=data[1:,6]
tank12_rna=array([tank12_rna_d[0],tank12_rna_d[2],tank12_rna_d[4],tank12_rna_d[6]])

tank3_rna_d=data[1:,7]
tank3_rna=array([tank3_rna_d[0],tank3_rna_d[2],tank3_rna_d[4],tank3_rna_d[6]])

tank8_rna_d=data[1:,8]
tank8_rna=array([tank8_rna_d[0],tank8_rna_d[2],tank8_rna_d[4],tank8_rna_d[6]])

tank10_rna_d=data[1:,9]
tank10_rna=array([tank10_rna_d[0],tank10_rna_d[2],tank10_rna_d[4],tank10_rna_d[6]])

tank4_rna_d=data[1:,10]
tank4_rna=array([tank4_rna_d[0],tank4_rna_d[2],tank4_rna_d[4],tank4_rna_d[6]])

tank7_rna_d=data[1:,11]
tank7_rna=array([tank7_rna_d[0],tank7_rna_d[2],tank7_rna_d[4],tank7_rna_d[6]])

tank11_rna_d=data[1:,12]
tank11_rna=array([tank11_rna_d[0],tank11_rna_d[2],tank11_rna_d[4],tank11_rna_d[6]])



data2=genfromtxt('Cell_density_final_eukaryotic.csv',delimiter=',') #cells/mL

day = data2[2:17,1]

### cell density ### now in cells/L

cd1_avv = 1000*data2[2:17,3]
cd1_av=array([cd1_avv[4],cd1_avv[7],cd1_avv[11],cd1_avv[14]])


cd2_avv = 1000*data2[2:17,7]
cd2_av=array([cd2_avv[4],cd2_avv[7],cd2_avv[11],cd2_avv[14]])

cd5_avv = 1000*data2[2:17,10]
cd5_av=array([cd5_avv[4],cd5_avv[7],cd5_avv[11],cd5_avv[14]])

cd6_avv =1000* data2[2:17,14]
cd6_av=array([cd6_avv[4],cd6_avv[7],cd6_avv[11],cd6_avv[14]])

cd9_avv = 1000*data2[2:17,18]
cd9_av=array([cd9_avv[4],cd9_avv[7],cd9_avv[11],cd9_avv[14]])

cd12_avv = 1000*data2[2:17,22]
cd12_av=array([cd12_avv[4],cd12_avv[7],cd12_avv[11],cd12_avv[14]])


data4=genfromtxt('cell_density_final_contaminant.csv',delimiter=',') 

day = data4[2:17,1]

### cell density ###
cd3_avv = 1000*data4[1:17,3]
cd3_av=array([cd3_avv[4],cd3_avv[7],cd3_avv[11],cd3_avv[14]])

cd4_avv = 1000*data4[1:17,7]
cd4_av=array([cd4_avv[4],cd4_avv[7],cd4_avv[11],cd4_avv[14]])

cd8_avv = 1000*data4[1:17,11]
cd8_av=array([cd8_avv[4],cd8_avv[7],cd8_avv[11],cd8_avv[14]])

cd7_avv = 1000*data4[1:17,15]
cd7_av=array([cd7_avv[4],cd7_avv[7],cd7_avv[11],cd7_avv[14]])

cd10_avv =1000* data4[1:17,19]
cd10_av=array([cd10_avv[4],cd10_avv[7],cd10_avv[11],cd10_avv[14]])

cd11_avv = 1000*data4[1:17,23]
cd11_av=array([cd11_avv[4],cd11_avv[7],cd11_avv[11],cd11_avv[14]])


Molar_mass_RNA_AT_average=261#316.47                        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=262#323.97   
Molar_mass_RNA=Molar_mass_RNA_AT_average*0.563+Molar_mass_RNA_CG_average*(1-0.563) 
Qc=1.00*10**(-12)/12

tank1_c_l=cd1_av*Qc
tank1_rna=(tank1_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank5_rna=(tank5_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank9_rna=(tank9_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank2_rna=(tank2_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank6_rna=(tank6_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank12_rna=(tank12_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
 
tank3_rna=(tank3_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank8_rna=(tank8_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank10_rna=(tank10_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank4_rna=(tank4_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank7_rna=(tank7_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)
tank11_rna=(tank11_rna)*(1e6)*(1/1e9)*(1/Molar_mass_RNA)#*(10.7/15.5)


CRNA1=cd1_avv*CRNA1*(1/10.7)
CRNA5=cd5_avv*CRNA5*(1/10.7)
CRNA9=cd9_avv*CRNA9*(1/10.7)
CRNA2=cd2_avv*CRNA2*(1/10.7)
CRNA6=cd6_avv*CRNA6*(1/10.7)
CRNA12=cd12_avv*CRNA12*(1/10.7)


CRNA3=cd3_avv[1:]*CRNA3*(1/10.7)
CRNA8=cd8_avv[1:]*CRNA8*(1/10.7)
CRNA10=cd10_avv[1:]*CRNA10*(1/10.7)
CRNA4=cd4_avv[1:]*CRNA4*(1/10.7)
CRNA7=cd7_avv[1:]*CRNA7*(1/10.7)
CRNA11=cd11_avv[1:]*CRNA11*(1/10.7)


days=arange(0,15,1)
day=[4,7,11,14]

figure(1,figsize=(8,6.5))

scatter(day,tank1_rna,color='#88CCEE', s=60, label='Tank 1')
scatter(day,tank5_rna,color='#228833', s=60, label='Tank 5')
scatter(day,tank9_rna,color='#4477AA', s=60, label='Tank 9')
scatter(day,tank2_rna,color='#DDCC77',marker='s', s=60, label='Tank 2')
scatter(day,tank6_rna,color='#CC6677',marker='s', s=60, label='Tank 6')
scatter(day,tank12_rna,color='#882255',marker='s', s=60, label='Tank 12')
scatter(day,tank3_rna,color='#88CCEE',marker='P', s=60, label='Tank 3')
scatter(day,tank8_rna,color='#228833',marker='P', s=60, label='Tank 8')
scatter(day,tank10_rna,color='#4477AA',marker='P', s=60, label='Tank 10')
scatter(day,tank4_rna,color='#DDCC77',marker='v', s=60, label='Tank 4')
scatter(day,tank7_rna,color='#CC6677',marker='v' ,s=60, label='Tank 7')
scatter(day,tank11_rna,color='#882255',marker='v' , s=60, label='Tank 11')
xlabel("Day of Experiment",fontsize=18)
title('Measured',fontsize=20)
ylabel('mol RNA L$^{-1}$',fontsize=18)


figure(2,figsize=(8,6.5))
scatter(days,CRNA1,color='#88CCEE', s=60, label='Tank 1')
scatter(days,CRNA5,color='#228833', s=60, label='Tank 5')
scatter(days,CRNA9,color='#4477AA', s=60, label='Tank 9')
scatter(days,CRNA2,color='#DDCC77',marker='s', s=60, label='Tank 2')
scatter(days,CRNA6,color='#CC6677',marker='s', s=60, label='Tank 6')
scatter(days,CRNA12,color='#882255',marker='s', s=60, label='Tank 12')
scatter(days,CRNA3,color='#88CCEE',marker='P', s=60, label='Tank 3')
scatter(days,CRNA8,color='#228833',marker='P', s=60, label='Tank 8')
scatter(days,CRNA10,color='#4477AA',marker='P', s=60, label='Tank 10')
scatter(days,CRNA4,color='#DDCC77',marker='v', s=60, label='Tank 4')
scatter(days,CRNA7,color='#CC6677',marker='v' ,s=60, label='Tank 7')
scatter(days,CRNA11,color='#882255',marker='v' , s=60, label='Tank 11')
title("Modeled",fontsize=20)
xlabel("Day of Experiment",fontsize=18)
ylabel(r'mol RNA L$^{-1}$',fontsize=18)


CRNA1=array([CRNA1[4],CRNA1[7],CRNA1[11],CRNA1[14]])
CRNA5=array([CRNA5[4],CRNA5[7],CRNA5[11],CRNA5[14]])
CRNA9=array([CRNA9[4],CRNA9[7],CRNA9[11],CRNA9[14]])
CRNA2=array([CRNA2[4],CRNA2[7],CRNA2[11],CRNA2[14]])
CRNA6=array([CRNA6[4],CRNA6[7],CRNA6[11],CRNA6[14]])
CRNA12=array([CRNA12[4],CRNA12[7],CRNA12[11],CRNA12[14]])


CRNA3=array([CRNA3[4],CRNA3[7],CRNA3[11],CRNA3[14]])
CRNA8=array([CRNA8[4],CRNA8[7],CRNA8[11],CRNA8[14]])
CRNA10=array([CRNA10[4],CRNA10[7],CRNA10[11],CRNA10[14]])
CRNA4=array([CRNA4[4],CRNA4[7],CRNA4[11],CRNA4[14]])
CRNA7=array([CRNA7[4],CRNA7[7],CRNA7[11],CRNA7[14]])
CRNA11=array([CRNA11[4],CRNA11[7],CRNA11[11],CRNA11[14]])

ex=arange(0,0.0005,0.00001)
ey=ex

figure(3,figsize=(8,6.5))
scatter(tank1_rna,CRNA1,color='#88CCEE', s=60, label='Tank 1',alpha=[0.25,0.5,0.75,1])
scatter(tank1_rna,CRNA1,facecolor='None', s=60, label='Tank 1',edgecolor='#88CCEE')
scatter(tank5_rna,CRNA5,color='#228833', s=60, label='Tank 5',alpha=[0.25,0.5,0.75,1])
scatter(tank5_rna,CRNA5,facecolor='None', s=60, label='Tank 5',edgecolor='#228833')
scatter(tank9_rna,CRNA9,color='#4477AA', s=60, label='Tank 9',alpha=[0.25,0.5,0.75,1])
scatter(tank9_rna,CRNA9,facecolor='None', s=60, label='Tank 9',edgecolor='#4477AA')
scatter(tank2_rna,CRNA2, color='#DDCC77', s=60,marker='s', label='Tank 2',alpha=[0.25,0.5,0.75,1])
scatter(tank2_rna,CRNA2, facecolor='None', s=60,marker='s', label='Tank 2',edgecolor='#DDCC77')
scatter(tank6_rna,CRNA6,color='#CC6677', s=60, marker='s',label='Tank 6',alpha=[0.25,0.5,0.75,1])
scatter(tank6_rna,CRNA6,facecolor='None', s=60,marker='s', label='Tank 6',edgecolor='#CC6677')
scatter(tank12_rna,CRNA12,color='#882255', s=60,marker='s', label= 'Tank 12',alpha=[0.25,0.5,0.75,1])
scatter(tank12_rna,CRNA12,facecolor='None', s=60,marker='s', label= 'Tank 12',edgecolor='#882255')
scatter(tank3_rna,CRNA3,color='#88CCEE', s=60,marker='P', label='Tank 3',alpha=[0,0.25,0.5,0.75,1])
scatter(tank3_rna,CRNA3,facecolor='None', s=60,marker='P', label='Tank 3',edgecolor='#88CCEE')
scatter(tank8_rna,CRNA8,color='#228833', s=60,marker='P', label='Tank 8',alpha=[0,0.25,0.5,0.75,1])
scatter(tank8_rna,CRNA8,facecolor='None', s=60,marker='P', label='Tank 8',edgecolor='#228833')
scatter(tank10_rna,CRNA10,color='#4477AA', s=60,marker='P', label='Tank 10',alpha=[0,0.25,0.5,0.75,1])
scatter(tank10_rna,CRNA10,facecolor='None', s=60,marker='P', label='Tank 10',edgecolor='#4477AA')
scatter(tank4_rna,CRNA4, color='#DDCC77', s=60,marker='v', label='Tank 4',alpha=[0,0.25,0.5,0.75,1])
scatter(tank4_rna,CRNA4, facecolor='None', s=60,marker='v', label='Tank 4',edgecolor='#DDCC77')
scatter(tank7_rna,CRNA7,color='#CC6677', s=60,marker='v', label='Tank 7',alpha=[0,0.25,0.5,0.75,1])
scatter(tank7_rna,CRNA7,facecolor='None', s=60,marker='v', label='Tank 7',edgecolor='#CC6677')
scatter(tank11_rna,CRNA11,color='#882255', s=60,marker='v', label= 'Tank 11',alpha=[0,0.25,0.5,0.75,1])
scatter(tank11_rna,CRNA11,facecolor='None', s=60,marker='v', label= 'Tank 11',edgecolor='#882255')
plot(ex,ey,color='black')
xlabel("Measured",fontsize=18)
ylabel('Modeled',fontsize=18)



show()
