
# PV=NRT
V = [2.25,5.50,1.40]
P=[250,250,760]
T=[20,20,8]
T1 = [0,0,0]
for j in range(len(T)):
    T1=T[j]+273
M = [2,32,28]
R = 0.082
m=[0,0,0]
for i in range(len(V)):
    m[i]=P[i]*V[i]*M[i]/R*T1[i]
print(m)
nomuser = input("donner le nom de l'element")
print(m[noms.indexe(nounuser)])


