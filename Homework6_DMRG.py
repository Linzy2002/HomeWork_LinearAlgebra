# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 20:56:15 2022

@author: lzy
"""
import numpy as np
import matplotlib.pyplot as plt
import math


def sortByEigenValue(Eigenvalues, EigenVectors):
#rearrange eigenvalues and eigenvectors by size of values
    index = np.argsort(Eigenvalues)
    Eigenvalues = Eigenvalues[index]
    EigenVectors = EigenVectors[:,index]
    return Eigenvalues, EigenVectors

#iteration means the max iterations,meanwhile It represents the number of nodes in the image.
Iteration=7

#retain express the matrix U and V retain size
retain=1000


#spin matrix
s1=np.matrix([[0,0.5],[0.5,0]])
s2=np.matrix([[0,-0.5*1j],[0.5*1j,0]])
sz=np.matrix([[0.5,0],[0,-0.5]])

ss=np.kron(s1,s1)+np.kron(s2,s2)+np.kron(sz,sz)

Hsu=np.kron(ss,np.eye(4))+np.kron(np.eye(2),np.kron(ss,np.eye(2)))+np.kron(np.eye(4),ss)


hs=ss
he=ss

e0=[0]

#The value when the loop approaches infinite length through iteration.
for k in range(1,Iteration+1):
    [v,e]=np.linalg.eig(Hsu)
    [vn,en]=sortByEigenValue(v,e)
    e0=np.append(e0,vn[0])

    m=np.sqrt(len(vn))
    m=int(m)
    rho=np.zeros((m,m))
    for i in range(1,m+1):
        rho[i-1,:]=en[(i-1)*m:i*m,0].T
    [U,S,V]=np.linalg.svd(rho)
    
    lenthU=len(S) 

#if the matrix is too big then retain the main component of it   
    if(lenthU>retain):
        U1=np.zeros((lenthU,retain))
        V1=np.zeros((lenthU,retain))
        for i in range(retain):
            U1[:,i]=U[:,i]
            V1[:,i]=V[:,i]
    else:
        U1=U
        V1=V
   
    h=len(hs)
    h05=int(h/2)
    
#update matrix for next iteration
    hs=np.kron(np.dot(np.dot(U1.T,hs),U1),np.eye(2)) + np.kron(np.dot(np.dot(U1.T,np.kron(np.eye(h05),s1)),U1),s1) + np.kron(np.dot(np.dot(U1.T,np.kron(np.eye(h05),s2)),U1),s2) + np.kron(np.dot(np.dot(U1.T,np.kron(np.eye(h05),sz)),U1),sz)
    he=np.kron(s1,np.dot(np.dot(V1.T,np.kron(s1,np.eye(h05))),V1)) + np.kron(s2,np.dot(np.dot(V1.T,np.kron(s2,np.eye(h05))),V1)) + np.kron(s2,np.dot(np.dot(V1.T,np.kron(s2,np.eye(h05))),V1)) + np.kron(np.eye(2),np.dot(np.dot(V1.T,he),V1))
    

    h2=len(hs)
    h205=int(h2/2)
    
    Hsu=np.kron(hs,np.eye(h2)) + np.kron(np.eye(h2),he) + np.kron(np.eye(h205),np.kron(ss,np.eye(h205)))
    print(e0)


lenthe=len(e0)

L=[]
eq=[]
st=[]

for i in range(lenthe):
    L=np.append(L,(i+1)*2+2)

for i in range(lenthe):
    eq=np.append(eq,e0[i]/L[i])
    st=np.append(st,0.25-math.log(2))
    
plt.plot(L,eq, 'ro-', alpha=0.8, linewidth=1, label='the answer is')
    
plt.plot(L,st)

plt.show()    
