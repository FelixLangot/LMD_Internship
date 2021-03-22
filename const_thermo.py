######################### .h file #############################
    #     ! const_thermo.h
    #     ! les principales constantes thermo utilisées dans calculs
        
    #     real CPD,CPV,CPVMCL,CL,RV,RD,LV0,G,ROWL,EPS,degK,
    #  :       Lfusion,lsub,lvap0


    #   parameter (CPD=1005.7) ! J/kg/K
    #   parameter (CPV=1870.0)      
    #   parameter (CL=4190.0)
    #   parameter (CPVMCL=CL-CPV)
    #   parameter (RV=461.5)
    #   parameter (RD=287.04)
    #   parameter (LV0=2.501E6)
    #   parameter (G=9.8)  
    #   parameter (ROWL=1000.0)
    #   parameter (EPS=RD/RV)
    #   parameter (degK=273.15)
    #   parameter (lsub=2.844E6) ! J/kg, from SAM
    #   parameter (lvap0=2.5104E6) ! J/kg, from SAM
    #   parameter (Lfusion=lsub - lvap0) ! J/kg, from SAM

##############################################################

Cpd = 1005.7 # J/kg/K      
Cpv = 1870
Cl = 4190
CpvMcl = Cl - Cpv
Rv = 461.5
Rd = 287.04
Lv0 = 2.501e6
g = 9.8
Rowl = 1000
Eps = Rd/Rv
degK = 273.15
lsub = 2.844e6       # J/kg, from SAM
lvap0 = 2.5104e6     # J/kg, from SAM
Lfusion = lsub - lvap0 # J/kg, from SAM
