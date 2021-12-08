import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
fE = open("E.dat","r")
E = fE.readlines()

nE = len(E)
for i in range(nE):
    E[i] = E[i][:-1].split(" ")
    E[i] = [int(s) for s in E[i] if s!=""]

fE.close()

fS = open("S.dat","r")
S = fS.readlines()

nS = len(S)
print(nS)
for i in range(nS):
    S[i] = S[i][:-1].split(" ")
    S[i] = [float(s) for s in S[i] if s!=""]


fig = plt.figure(figsize=(10,10)) 
ax = Axes3D(fig)
ax.set_axis_off()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)


Vtot=0

for i in range(nE):
    triangle = E[i]
    x = [S[sommet][0] for sommet in triangle]
    y = [S[sommet][1] for sommet in triangle]
    z = [S[sommet][2] for sommet in triangle]
    
    A=np.array([x[0],y[0],z[0]])
    B=np.array([x[1],y[1],z[1]])
    C=np.array([x[2],y[2],z[2]])
    D=np.array([0,0,0])

    A=A.reshape(3,1)
    B=B.reshape(3,1)
    C=C.reshape(3,1)
    D=D.reshape(3,1)

    AB=B-A
    AC=C-A
    AD=D-A
    V=np.zeros((3,3))
    for i in range(3):
        V[i,0]=AB[i]
        V[i,1]=AC[i]
        V[i,2]=AD[i]
    V=1/6*abs(np.linalg.det(V))
    Vtot+=V

    verts = [list(zip(x,y,z))]
    tri = Poly3DCollection(verts)
    tri.set_color([1,0.5,0.5])
    tri.set_edgecolor('k')
    
    ax.add_collection3d(tri)
plt.savefig('sphere_'+str(nS)+'.png')
#plt.show()

print("Volumee approx: ",Vtot)
print("Vrai volume: ",4/3*3.14159)
print("Erreur: ",4/3*3.14159-Vtot)
