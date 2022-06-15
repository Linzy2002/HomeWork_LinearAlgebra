# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 09:38:03 2022

@author: lzy
"""

import numpy as np
import matplotlib.pyplot as plt


def producef0(D):
    f0=np.random.rand(D)
    f0=f0/np.linalg.norm(f0)
    return f0

def produceH(L,delta):
    m=2**L
    
    s1=np.matrix([[0,1],[0,0]])
    s2=np.matrix([[0,0],[1,0]])
    sz=np.matrix([[0.5,0],[0,-0.5]])
    
    s3=np.kron(s1,s2)
    s4=np.kron(s2,s1)
    sz2=np.kron(sz,sz)
    hi=0.5*(s3+s4)+delta*sz2
    H=np.zeros([m,m])
    
    for i in range(L-1):
        H+=np.kron(np.kron(np.eye(2**i),hi),np.eye(2**(L-i-2)))
    
    H+=0.5*np.kron(np.kron(s2,np.eye(2**(L-2))),s1)  
    H+=0.5*np.kron(np.kron(s1,np.eye(2**(L-2))),s2)
    H+=0.5*np.kron(np.kron(sz,np.eye(2**(L-2))),sz)
    
    return H


def lanczos(H,f0,f_1,i,HD,HS):
    Hf0=np.dot(H,f0)
    hd =np.dot(f0,Hf0)
    h10=np.dot(f_1,Hf0)
    f1 = Hf0 - hd * f0 - h10 * f_1
    f1=f1/np.linalg.norm(f1)
    hs = np.dot(f1,Hf0)
    HD=np.append(HD,hd)
    HS=np.append(HS,hs)
    #Hf=makeHf(i+1,HD,HS)
    
    Hf=np.zeros([i+1,i+1])
    for j in range(i+1):
        Hf[j,j]=HD[j]
    for j in range(i):
        Hf[j+1,j]=HS[j]
        Hf[j,j+1]=HS[j]
    
    e,v=np.linalg.eig(Hf)
    index=np.argmin(e)
    E0=e[index]
    
    return f1,E0,HD,HS


    
def SolveEnergy(L,delta,err):
    m=2**L
    Eg=0

    f_1=np.zeros(m)
    f0=producef0(m)
    
    HD=[]
    HS=[]
    H=produceH(L,delta)
    for i in range(m):
        [f1,E0,HDn,HSn]=lanczos(H,f0,f_1,i,HD,HS)
        
        if(abs(E0-Eg)<err):
            break
        
        Eg=E0
        f_1=f0
        f0=f1
        HD=HDn
        HS=HSn
    return Eg

def FindPoint(L,d1,d2,T):
    err=0.0001
    
    for i in range(T):
        energy=[]
        
        delta=np.linspace(d1,d2,11)
    
        for d in delta:
            energy=np.append(energy,SolveEnergy(L,d,err))

        d_ans=d1+0.1*np.argmax(energy)*(d2-d1)
        
        dr=d2-d1
        
        d1=d_ans-0.1*dr
        d2=d_ans+0.1*dr
        print("now is finding phase trasition point at precision",d2-d1,"\n")
        
        precision=d2-d1
        
    return d_ans,precision



        

if __name__=="__main__":
    err=0.0001
    L=12
    energy=[]
    delta=list(np.arange(-1.5,2,0.5))
    for d in delta:
        energy=np.append(energy,SolveEnergy(L,d,err))
        print("now is computing energy at delta =",d)

    plt.plot(delta, energy, 'ro-', alpha=0.8, linewidth=1, label='the answer is')

    T=3
    
    E,precision=FindPoint(L,-1.5,-0.5,T)
    print("The phase transition point is at the delta of ",E,"at the precision of ",precision)
