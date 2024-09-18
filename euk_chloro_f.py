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




data2=genfromtxt('Cell_density_final_eukaryotic.csv',delimiter=',') #cells/mL

day = data2[2:17,1]

### cell density ### now in cells/L

cd1_av = 1000*data2[2:17,3]
cd1_av=array([cd1_av[4],cd1_av[7],cd1_av[11],cd1_av[14]])

cd1_upper = 1000*data2[2:17,4]
cd1_lower = 1000*data2[2:17,5]

cd2_av = 1000*data2[2:17,7]
cd2_av=array([cd2_av[4],cd2_av[7],cd2_av[11],cd2_av[14]])

cd2_upper = 1000*data2[2:17,8]
cd2_lower = 1000*data2[2:17,9]

cd5_av = 1000*data2[2:17,10]
cd5_av=array([cd5_av[4],cd5_av[7],cd5_av[11],cd5_av[14]])

cd5_upper = 1000*data2[2:17,12]
cd5_lower = 1000*data2[2:17,13]

cd6_av =1000* data2[2:17,14]
cd6_av=array([cd6_av[4],cd6_av[7],cd6_av[11],cd6_av[14]])

cd6_upper = 1000*data2[2:17,16]
cd6_lower =1000* data2[2:17,17]

cd9_av = 1000*data2[2:17,18]
cd9_av=array([cd9_av[4],cd9_av[7],cd9_av[11],cd9_av[14]])

cd9_upper =1000* data2[2:17,20]
cd9_lower =1000* data2[2:17,21]

cd12_av = 1000*data2[2:17,22]
cd12_av=array([cd12_av[4],cd12_av[7],cd12_av[11],cd12_av[14]])

cd12_upper =1000* data2[2:17,24]
cd12_lower =1000* data2[2:17,25]

data3 = genfromtxt('chlorophyll.csv',delimiter=',')

T1_chl = data3[2:17,1]
T2_chl = data3[2:17,3]
T5_chl = data3[2:17,5]
T6_chl = data3[2:17,7]
T6_chl=T6_chl[0:3]
T9_chl = data3[2:17,9]
T12_chl = data3[2:17,11]

mT1_chl = (59/55)*(893.51)*(1e6)*cd1_av*array([Cchl1[4],Cchl1[7],Cchl1[11],Cchl1[14]])
mT2_chl = (59/55)*(893.51)*(1e6)*cd2_av*array([Cchl2[4],Cchl2[7],Cchl2[11],Cchl2[14]])
mT5_chl = (59/55)*(893.51)*(1e6)*cd5_av*array([Cchl5[4],Cchl5[7],Cchl5[11],Cchl5[14]])
mT6_chl = (59/55)*(893.51)*(1e6)*cd6_av*array([Cchl6[4],Cchl6[7],Cchl6[11],Cchl6[14]])
mT6_chl=mT6_chl[0:3]
mT9_chl = (59/55)*(893.51)*(1e6)*cd9_av*array([Cchl9[4],Cchl9[7],Cchl9[11],Cchl9[14]])
mT12_chl = (59/55)*(893.51)*(1e6)*cd12_av*array([Cchl12[4],Cchl12[7],Cchl12[11],Cchl12[14]])

data4=genfromtxt('cell_density_final_contaminant.csv',delimiter=',') 

day = data4[2:17,1]

### cell density ###

cd3_av = 1000*data4[1:17,3]
cd3_av=array([cd3_av[4],cd3_av[7],cd3_av[11],cd3_av[14]])



cd4_av = 1000*data4[1:17,7]
cd4_av=array([cd4_av[4],cd4_av[7],cd4_av[11],cd4_av[14]])

cd8_av = 1000*data4[1:17,11]
cd8_av=array([cd8_av[4],cd8_av[7],cd8_av[11],cd8_av[14]])

cd7_av = 1000*data4[1:17,15]
cd7_av=array([cd7_av[4],cd7_av[7],cd7_av[11],cd7_av[14]])

cd10_av =1000* data4[1:17,19]
cd10_av=array([cd10_av[4],cd10_av[7],cd10_av[11],cd10_av[14]])

cd11_av = 1000*data4[1:17,23]
cd11_av=array([cd11_av[4],cd11_av[7],cd11_av[11],cd11_av[14]])

mT3_chl = (59/55)*(893.51)*(1e6)*cd3_av*array([Cchl3[4],Cchl3[7],Cchl3[11],Cchl3[14]])
mT4_chl = (59/55)*(893.51)*(1e6)*cd4_av*array([Cchl4[4],Cchl4[7],Cchl4[11],Cchl4[14]])
mT8_chl = (59/55)*(893.51)*(1e6)*cd8_av*array([Cchl8[4],Cchl8[7],Cchl8[11],Cchl8[14]])
mT7_chl = (59/55)*(893.51)*(1e6)*cd7_av*array([Cchl7[4],Cchl7[7],Cchl7[11],Cchl7[14]])
mT10_chl = (59/55)*(893.51)*(1e6)*cd10_av*array([Cchl10[4],Cchl10[7],Cchl10[11],Cchl10[14]])
mT11_chl = (59/55)*(893.51)*(1e6)*cd11_av*array([Cchl11[4],Cchl11[7],Cchl11[11],Cchl11[14]])

data5 = genfromtxt('chlorophyll_final.csv',delimiter=',')

T3_chl = data5[7,3:7]
T8_chl = data5[8,3:7]
T10_chl = data5[9,3:7]
T4_chl = data5[10,3:7]
T7_chl = data5[11,3:7]
T11_chl = data5[12,3:7]

mT4_chl=delete(mT4_chl,2)
T4_chl=delete(T4_chl,2)

x=append(T1_chl,T5_chl)
x2=append(T9_chl,T3_chl)
#x3=append(x,x2)
# x4=append(T6_chl,T12_chl)
# x5=append(T3_chl,T4_chl)
x6=append(x,x2)
x7=append(T8_chl,T10_chl)
# x8=append(T10_chl,T11_chl)
x9=append(x6,x7)
#x10=append(x3,x6)
x11=x9#append(x10,x9)
  
y=append(mT1_chl,mT5_chl)
y2=append(mT9_chl,mT3_chl)
y3=append(y,y2)
# y4=append(mT6_chl,mT12_chl)
# y5=append(mT3_chl,mT4_chl)
# y6=append(y4,y5)
y7=append(mT8_chl,mT10_chl)
# y8=append(mT10_chl,mT11_chl)
# y9=append(y7,y8)
y10=append(y3,y7)
y11=y10#append(y10,y9)



y11.sort()
x11.sort()
model=LinearRegression()
  
model_fit=model.fit(x11.reshape(-1,1),y11)
rq_sq=model.score(x11.reshape(-1,1),y11)
y_pred=model.predict(x11.reshape(-1,1))



reg_leg=mpat.Patch(color='#fff',label='Regression Line:',alpha=0.75)
slope_leg=mpat.Patch(color='#fff',label='Slope=1.67',alpha=0.75)
r2_leg=mpat.Patch(color='#fff',label='R$^{2}$=0.88',alpha=0.75)
handles=[reg_leg,slope_leg,r2_leg]



ex=arange(0,351,1)
ey=ex

figure(1,figsize=(8,6.5))
scatter(T1_chl,mT1_chl,color='#88CCEE', s=60, label='Tank 1',alpha=[0.25,0.5,0.75,1])
scatter(T1_chl,mT1_chl,facecolor='None', s=60, label='Tank 1',edgecolor='#88CCEE')
scatter(T5_chl,mT5_chl,color='#228833', s=60, label='Tank 5',alpha=[0.25,0.5,0.75,1])
scatter(T5_chl,mT5_chl,facecolor='None', s=60, label='Tank 5',edgecolor='#228833')
scatter(T9_chl,mT9_chl,color='#4477AA', s=60, label='Tank 9',alpha=[0.25,0.5,0.75,1])
scatter(T9_chl,mT9_chl,facecolor='None', s=60, label='Tank 9',edgecolor='#4477AA')
scatter(T3_chl,mT3_chl,color='#88CCEE', s=60,marker='P', label='Tank 3',alpha=[0,0.25,0.5,0.75,1])
scatter(T3_chl,mT3_chl,facecolor='None', s=60,marker='P', label='Tank 3',edgecolor='#88CCEE')
scatter(T8_chl,mT8_chl,color='#228833', s=60,marker='P', label='Tank 8',alpha=[0,0.25,0.5,0.75,1])
scatter(T8_chl,mT8_chl,facecolor='None', s=60,marker='P', label='Tank 8',edgecolor='#228833')
scatter(T10_chl,mT10_chl,color='#4477AA', s=60,marker='P', label='Tank 10',alpha=[0,0.25,0.5,0.75,1])
scatter(T10_chl,mT10_chl,facecolor='None', s=60,marker='P', label='Tank 10',edgecolor='#4477AA')
plot(x11,y_pred,color='black',linestyle='dashed')

plot(ex,ey,color='black')
xlabel('Measured',fontsize=18)
ylabel('Modeled',fontsize=18)
 
xt=arange(0,400,50)
xticks(ticks=xt,fontsize=16)
yticks(ticks=xt,fontsize=16)
ylim(-5,350)
xlim(-5,350)
legend(handles=handles,loc='lower right',fontsize='x-large',frameon=False)
title('Chlorophyll-$\it{a}$ ($\mu$g L$^{-1}$)',fontsize=20)




x=append(T1_chl,T5_chl)
x2=append(T9_chl,T2_chl)
x3=append(x,x2)
x4=append(T6_chl,T12_chl)
x5=append(T3_chl,T4_chl)
x6=append(x4,x5)
x7=append(T8_chl,T7_chl)
x8=append(T10_chl,T11_chl)
x9=append(x7,x8)
x10=append(x3,x6)
x11=append(x10,x9)
   
y=append(mT1_chl,mT5_chl)
y2=append(mT9_chl,mT2_chl)
y3=append(y,y2)
y4=append(mT6_chl,mT12_chl)
y5=append(mT3_chl,mT4_chl)
y6=append(y4,y5)
y7=append(mT8_chl,mT7_chl)
y8=append(mT10_chl,mT11_chl)
y9=append(y7,y8)
y10=append(y3,y6)
y11=append(y10,y9)

y11.sort()
x11.sort()
model=LinearRegression()
  
model_fit=model.fit(x11.reshape(-1,1),y11)
rq_sq=model.score(x11.reshape(-1,1),y11)
y_pred=model.predict(x11.reshape(-1,1))


   
print(rq_sq)
print(f"intercept: {model.intercept_}")
print(f"slope: {model.coef_}")

reg_leg=mpat.Patch(color='#fff',label='Regression Line:',alpha=0.75)
slope_leg=mpat.Patch(color='#fff',label='Slope=1.87',alpha=0.75)
r2_leg=mpat.Patch(color='#fff',label='R$^{2}$=0.93',alpha=0.75)
handles=[reg_leg,slope_leg,r2_leg]



ex=arange(0,351,1)
ey=ex

figure(2,figsize=(8,6.5))
scatter(T1_chl,mT1_chl,color='#88CCEE', s=60, label='Tank 1',alpha=[0.25,0.5,0.75,1])
scatter(T1_chl,mT1_chl,facecolor='None', s=60, label='Tank 1',edgecolor='#88CCEE')
scatter(T5_chl,mT5_chl,color='#228833', s=60, label='Tank 5',alpha=[0.25,0.5,0.75,1])
scatter(T5_chl,mT5_chl,facecolor='None', s=60, label='Tank 5',edgecolor='#228833')
scatter(T9_chl,mT9_chl,color='#4477AA', s=60, label='Tank 9',alpha=[0.25,0.5,0.75,1])
scatter(T9_chl,mT9_chl,facecolor='None', s=60, label='Tank 9',edgecolor='#4477AA')
scatter(T2_chl,mT2_chl, color='#DDCC77', s=60,marker='s', label='Tank 2',alpha=[0.25,0.5,0.75,1])
scatter(T2_chl,mT2_chl, facecolor='None', s=60,marker='s', label='Tank 2',edgecolor='#DDCC77')
scatter(T6_chl,mT6_chl,color='#CC6677', s=60, marker='s',label='Tank 6',alpha=[0.25,0.5,0.75,1])
scatter(T6_chl,mT6_chl,facecolor='None', s=60,marker='s', label='Tank 6',edgecolor='#CC6677')
scatter(T12_chl,mT12_chl,color='#882255', s=60,marker='s', label= 'Tank 12',alpha=[0.25,0.5,0.75,1])
scatter(T12_chl,mT12_chl,facecolor='None', s=60,marker='s', label= 'Tank 12',edgecolor='#882255')
scatter(T3_chl,mT3_chl,color='#88CCEE', s=60,marker='P', label='Tank 3',alpha=[0,0.25,0.5,0.75,1])
scatter(T3_chl,mT3_chl,facecolor='None', s=60,marker='P', label='Tank 3',edgecolor='#88CCEE')
scatter(T8_chl,mT8_chl,color='#228833', s=60,marker='P', label='Tank 8',alpha=[0,0.25,0.5,0.75,1])
scatter(T8_chl,mT8_chl,facecolor='None', s=60,marker='P', label='Tank 8',edgecolor='#228833')
scatter(T10_chl,mT10_chl,color='#4477AA', s=60,marker='P', label='Tank 10',alpha=[0,0.25,0.5,0.75,1])
scatter(T10_chl,mT10_chl,facecolor='None', s=60,marker='P', label='Tank 10',edgecolor='#4477AA')
scatter(T4_chl,mT4_chl, color='#DDCC77', s=60,marker='v', label='Tank 4',alpha=[0,0.25,0.5,0.75,1])
scatter(T4_chl,mT4_chl, facecolor='None', s=60,marker='v', label='Tank 4',edgecolor='#DDCC77')
scatter(T7_chl,mT7_chl,color='#CC6677', s=60,marker='v', label='Tank 7',alpha=[0,0.25,0.5,0.75,1])
scatter(T7_chl,mT7_chl,facecolor='None', s=60,marker='v', label='Tank 7',edgecolor='#CC6677')
scatter(T11_chl,mT11_chl,color='#882255', s=60,marker='v', label= 'Tank 11',alpha=[0,0.25,0.5,0.75,1])
scatter(T11_chl,mT11_chl,facecolor='None', s=60,marker='v', label= 'Tank 11',edgecolor='#882255')
plot(x11,y_pred,color='black',linestyle='dashed')

plot(ex,ey,color='black')
xlabel('Measured',fontsize=18)
ylabel('Modeled',fontsize=18)
 
xt=arange(0,400,50)
xticks(ticks=xt,fontsize=16)
yticks(ticks=xt,fontsize=16)
ylim(-5,350)
xlim(-5,350)
legend(handles=handles,loc='lower right',fontsize='x-large',frameon=False)
title('Chlorophyll-$\it{a}$ ($\mu$g L$^{-1}$)',fontsize=20)




show()
