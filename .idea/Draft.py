import chemicals
from thermo.chemical import Chemical
from thermo import VaporPressure
from scipy.optimize import fsolve
from sympy import *
import numpy as np
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
n_Pentane= Chemical('n-Pentane')
print(n_Pentane.Tc)
print(n_Pentane.Tb)
print(n_Pentane.Pc)
n_Pentane_psat = VaporPressure(Tb=309.21 , Tc=469.7, Pc= 3370000.0 , CASRN='109-66-0')
print(n_Pentane_psat.all_methods)
b=chemicals.vapor_pressure.Psat_data_AntoinePoling
print(b)
print(chemicals.vapor_pressure.Psat_data_AntoinePoling.axes)
#PRINT THE HEAD OF THE DataFrame
print(chemicals.vapor_pressure.Psat_data_AntoinePoling.head(5))
#PRINT THE COLUMN OF DATAFRAME
print(chemicals.vapor_pressure.Psat_data_AntoinePoling ['A'])
#PRINT LINGE OF DATFRAME
print(chemicals.vapor_pressure.Psat_data_AntoinePoling.iloc[121])
n_Butane= Chemical('n-Butane')
print(n_Butane.Tb)
print(n_Butane.Tc)
print(n_Butane.Pc)
n_Butane_psat = VaporPressure(Tb= 272.65, Tc=425.12, Pc= 3796000.0 , CASRN='106-97-8')
print(chemicals.vapor_pressure.Psat_data_AntoinePoling.iloc[92])
Propane= Chemical('Propane')
print(Propane.Tb)
print(Propane.Tc)
print(Propane.Pc)
Propane_psat = VaporPressure(Tb= 231.04, Tc=369.83, Pc= 4248000.0 , CASRN='74-98-6')
print(chemicals.vapor_pressure.Psat_data_AntoinePoling.iloc[21])
print(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21])
print(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])
print(chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])

x=[[0.8894143384480268, 0.09713206869323096, 0.013453592858742194], [0.6258019950213685, 0.2541628135964542, 0.12003519138217728], [0.34443944979040364, 0.2311262047713737, 0.42443434543822256], [0.24274673849284256, 0.16288823020307522, 0.5943650313040821], [0.2427467384928426, 0.16288823020307522, 0.5943650313040822]]
#Température de l'étage 1
def func(T):
   for i in range(1):
    X= 10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])+T))/(689476)*x[0][i]
    Y=10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[92]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[92])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[92])+T))/(689476)*x[0][i]
    Z=10**( chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[121]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[121])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[121])+T))/(689476)*x[0][i]
    S=X+Y+Z-1
   return S
root=fsolve(func,291)
print("T étage1:",root)
#Température de l'étage 2
def func(T):
   for i in range(1):
    X= 10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])+T))/(689476)*x[1][i]
    Y=10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[92]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[92])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[92])+T))/(689476)*x[1][i]
    Z=10**( chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[121]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[121])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[121])+T))/(689476)*x[1][i]
    S=X+Y+Z-1
   return S
root1=fsolve(func,291)
print("T étage2:",root1)
#Température de l'étage 3
def func(T):
   for i in range(1):
    X= 10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])+T))/(689476)*x[2][i]
    Y=10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[92]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[92])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[92])+T))/(689476)*x[2][i]
    Z=10**( chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[121]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[121])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[121])+T))/(689476)*x[2][i]
    S=X+Y+Z-1
   return S
root2=fsolve(func,291)
print("T étage3:",root2)
#Température de l'étage 4
def func(T):
   for i in range(1):
    X= 10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])+T))/(689476)*x[3][i]
    Y=10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[92]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[92])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[92])+T))/(689476)*x[3][i]
    Z=10**( chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[121]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[121])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[121])+T))/(689476)*x[3][i]
    S=X+Y+Z-1
   return S
root3=fsolve(func,291)
print("T étage4:",root3)
#Température de l'étage 5
def func(T):
   for i in range(1):
    X= 10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[21]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[21])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[21])+T))/(689476)*x[4][i]
    Y=10**(chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[92]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[92])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[92])+T))/(689476)*x[4][i]
    Z=10**( chemicals.vapor_pressure.Psat_data_AntoinePoling['A'].iloc[121]-(chemicals.vapor_pressure.Psat_data_AntoinePoling['B'].iloc[121])/((chemicals.vapor_pressure.Psat_data_AntoinePoling['C'].iloc[121])+T))/(689476)*x[4][i]
    S=X+Y+Z-1
   return S
root4=fsolve(func,291)
print("T étage5:",root4)











