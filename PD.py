#Head of code
import numpy as np
from thermo.chemical import Chemical
from thermo import VaporPressure
from scipy.optimize import fsolve
import chemicals
#input
T= [291, 305, 319, 333, 347]
P=np.array([689476,689476,689476,689476,689476])
p=[0,0,0,0,0]
F=[0,0,100,0,0]
U=[50,0,0,0,0]
V= np.array([0,150,150,150,150,0])
D= np.array([0,0,-30,0,0])
#CODE POUR PROPANE
#identification de l'élément
Propane_psat = VaporPressure(Tb= 231.04, Tc=369.83, Pc= 4248000.0 , CASRN='74-98-6')
#calcule de pression vapeur saturante pour chaque T
for i in range(len(T)):
   p[i]=Propane_psat.calculate(T=T[i], method='ANTOINE_POLING')
#Calcule de volatilité
K3=[]
for i in range(5):
    M= p[i]/P[i]
    K3.append(M)
print("les valeur de la volatilité propane : ")
print(K3)
#calcule de la matrice
A=[0,0,0,0,0]
B=[0,0,0,0,0]
C=[0,0,0,0,0]
S1= [0,0,0,0,0]
#CALCULE DES SOMMES
for j in range(5):
    S1[j]=S1[j-1]+(F[j-1]-U[j-1])
S2= [0,0,0,0,0]
for j in range(4):
    S2[j]=S2[j-1]+(F[j]-U[j])
#calcule de la  matrice
for i in range(5):
    B[i]=-V[i+1]-S2[i]-U[i]-V[i]*K3[i]
    A[i]=V[i]+S1[i]
for i in range(4):
    C[i]= V[i+1]*K3[i+1]
C3= [[B[0],C[0],0,0,0],
 [A[1],B[1],C[1],0,0],
 [0,A[2],B[2],C[2],0],
 [0,0,A[3],B[3],C[3]],
 [0,0,0,A[4],B[4]]]
print("la matrice de Propane" )
print(C3)
#résolution
X = np.linalg.inv(C3).dot(D)
print("la valeur de x : ")
print(X)
#code pur n_butane
n_Butane_psat = VaporPressure(Tb= 272.65, Tc=425.12, Pc= 3796000.0 , CASRN='106-97-8')
for i in range(len(T)):
    p[i]=n_Butane_psat.calculate(T=T[i],method='ANTOINE_POLING')
#calcule de volatilité
K4=[]
for i in range(5):
    M= p[i]/P[i]
    K4.append(M)
print("les valeur de la volatilité de n_Butane : ")
print(K4)
#calcule de la matrice
for i in range(1):
    B[i]=-V[i+1]-S2[i]-U[i]-V[i]*K4[i]
    A[i]=V[i]+S1[i]
    C[i]= V[i+1]*K4[i+1]
C4= [[B[0],C[0],0,0,0],
 [A[1],B[1],C[1],0,0],
 [0,A[2],B[2],C[2],0],
 [0,0,A[3],B[3],C[3]],
 [0,0,0,A[4],B[4]]]
print("la matrice de n_Butane ")
print(C4)
X4 = np.linalg.inv(C4).dot(D)
print("la valeur de x : ")
print(X4)
#code pur n_Pentane
n_Pentane_psat = VaporPressure(Tb=309.21 , Tc=469.7, Pc= 3370000.0 , CASRN='109-66-0')
for i in range(len(T)):
    p[i]=n_Pentane_psat.calculate(T=T[i],method='ANTOINE_POLING')
P=np.array([689476,689476,689476,689476,689476])
K5=[]
for i in range(5):
    M= p[i]/P[i]
    K5.append(M)
print("les valeur de la volatilité de n_Pentane  : ")
print(K5)
#calcule de la matrice
for i in range(4):
    B[i]=-V[i+1]-S2[i]-U[i]-V[i]*K5[i]
    A[i]=V[i]+S1[i]
    C[i]= V[i+1]*K5[i+1]
C5= [[B[0],C[0],0,0,0],
 [A[1],B[1],C[1],0,0],
 [0,A[2],B[2],C[2],0],
 [0,0,A[3],B[3],C[3]],
 [0,0,0,A[4],B[4]]]
print("la matrice de n_Pentane ")
print(C5)
D5=np.array([0,0,-40,0,0])
X5 = np.linalg.inv(C5).dot(D5)
print("la valeur de x : ")
print(X5)
SX=[]
for i in range(5):
    a=X[i]+X4[i]+X5[i]
    SX.append(a)
print("la somme des élément pour chaque étage ",SX)
#normalisation
E1=[0,0,0,0,0]
E2=[0,0,0,0,0]
E3=[0,0,0,0,0]
for j in range(5):
    E1[j]= X[j]/SX[j]
    E2[j]= X4[j]/SX[j]
    E3[j]=X5[j]/SX[j]
x= [[E1[0],E2[0],E3[0]],
 [E1[1],E2[1],E3[1]],
 [E1[2],E2[2],E3[2]],
 [E1[3],E2[3],E3[3]],
 [E1[4],E2[4],E3[4]]]
print("la matrice des éléments pour chaque étage")
print(x)
somme=[0,0,0,0,0]
for i in range(5):
  somme[i]=E1[i]+E2[i]+E3[i]
print("la somme des élément pour chaque étage apres normalisation",somme)
#construction de matrice K de meme dimension que x
K=[0,0,0,0,0]
for i in range(5):
  K[i]=[K3[i],K4[i],K5[i]]
print(K)
sum_étage1= np.dot(K[0],x[0])
print(sum_étage1)
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
#nouvaux vecteur de température
T=np.array([root,root1,root2,root3,root4])
print(np.array(T))






