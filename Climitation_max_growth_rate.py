'''
Created on May 8, 2024

@author: garmin
'''
from numpy import *
from af001_energy_calculation import *

#CFM constants either from Kei's code emailed to me or supplemental table in inomura 2020
def Clim_mu(I):
#light and photosynthetic constants/equations
    I=I
    E3=evalue()
    #E=E3.E
    E=0.7741553538213819
    
    #m=3.93e-1 #3.79146798299876E-19                                #3.79146798299876E-19
    m=0.3930994004773114
    vIm=276.923478810005 #2.77e2 #0.00320513285659728                          #0.00320513285659728        #maximum photosynthetic rate called Pmax in OG CFM code
    AI= 0.00863364097132997                           #0.00863364097132997         #called OT in CFM
    vI=vIm*(1-exp(-AI*I))
    vI=vIm 
    Achl=(1+E)/vI
    Bchl=m/vI
     
    #proportionality constants
    #Apho= 1.60e1 # 3.56099164557551          #1.60e1  
    Apho=15.988852488634041     
    #Abio=  2.71e-1 # 4.34728279914354E-10   # 2.71e-1  
    Abio=0.2711013856688124
    #AP_RNA= 4.23e-3#6212.59249917364     # 4.23e-3
    AP_RNA=0.004234953828227311
    A_pchl_pho= 0.0281633095303638   
     
    #stoichiometric constants
    CP_Plip=40
    CP_RNA=10.7
    NP_RNA=3.8
    NC_Chl=4/55
    NC_pro=1/4.49
    NC_DNA=3.8/11.1
    
    
    
        
    Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
    
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    
    YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    
    CP_RNA=YrnaC_N/YnucacidP_N
    
    
     
    #macromolecular constants
    Qp_RNA_min=0.00022282635312159368#2.23e-4             #got this from inomura 2020 supplemental
    Qc_DNA=0.0009414247529173031#9.41e-4                            #got this from inomura 2020 supplemental
    
    Nconst_protein=4.45336898828389E-15 
    Qc=1.00*10**(-12)/12 
    Qc_oth=0.01821441041892576#1.51786753491048E-15/Qc     #Cessential in CFM
    Qc_pro_oth=0.239947521088736#Nconst_protein/Qc*(1/NC_pro)
    Qc_other=Qc_pro_oth+Qc_DNA+Qc_oth
    
    #Climitation
    a=CP_RNA*AP_RNA*((Apho*Achl)+Abio)
    b=Achl*(1+Apho+CP_Plip*A_pchl_pho)+CP_RNA*AP_RNA*(Apho*Bchl+Qc_pro_oth)+Abio
    c=Bchl*(1+Apho+CP_Plip*A_pchl_pho)+CP_RNA*Qp_RNA_min+Qc_other-1
    
    #Nlimitation
    #a=NP_RNA*AP_RNA*((Apho*Achl)+Abio)
    #b=NC_Chl*Achl+NC_pro*(Apho*Achl+Abio)+NP_RNA*AP_RNA*(Apho*Bchl+Qc_pro_oth)
    #c=NC_Chl*Bchl+NC_pro*(Apho*Bchl+Qc_pro_oth)+NC_DNA*Qc_DNA+NP_RNA*Qp_RNA_min-Qn
    
    mu_max=(-b+(b**2-(4*a*c))**(1/2))/(2*a)
    
    return(mu_max)













