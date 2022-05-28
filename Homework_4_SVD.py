# -*- coding: utf-8 -*-
"""
Created on Sat May 28 18:05:37 2022

@author: lzy
"""

import numpy as np

def produce(m,n):
    a=np.random.rand(m,n)
    b=np.random.rand(m,n)
    a=a+1j*b
    return a

def sortByEigenValue(Eigenvalues, EigenVectors):
#rearrange eigenvalues and eigenvectors by size of values
    index = np.argsort(-1*Eigenvalues)
    Eigenvalues = Eigenvalues[index]
    EigenVectors = EigenVectors[:,index]
    return Eigenvalues, EigenVectors


def SVD( a ):
    (m, n) = np.shape(a)  # get the size of matrix a 
    r = np.linalg.matrix_rank(a) #get the number of singular in matrix a
    
    
    if (m<n):        #to symply the method of caculate,judge the shape 
        
        adaggera=a.conjugate().T.dot(a)
    
        value, v = np.linalg.eig(adaggera)
        value2,v2=sortByEigenValue(value, v) 
        #rearrange eigenvalues and eigenvectors
        
        s = np.sqrt(value2)[:r].real
        u = np.zeros((m,m),dtype='complex')
    
    
        for i in range(r):
            u[:,i] = np.dot(a,v2[:,i])/s[i]
        for i in range(r,m):
            v1 = np.random.rand(m)
            for j in range(i):
                v1 -= np.dot(u[:,j].conjugate(),v1) * u[:,j]
                v1 /= np.linalg.norm(v1)
                u[:,i] = v1
        
        S=np.zeros([m,n])
    
        for i in range(r):
            S[i,i]=s[i]
        
        return u, S, v2.conjugate().T


    
    else :
        
        aadagger=a.dot(a.conjugate().T)
        
        value, v = np.linalg.eig(aadagger)
        
        value2,v2=sortByEigenValue(value, v)
        s = np.sqrt(value2)[:r].real
        u = np.zeros((n,n),dtype='complex')
        
        for i in range(r):
            u[:,i] = np.dot(a.conjugate().T,v2[:,i])/s[i]
        for i in range(r,n):
            v1 = np.random.rand(n)
            for j in range(i):
                v1 -= np.dot(u[:,j].conjugate(),v1) * u[:,j]
                v1 /= np.linalg.norm(v1)
                u[:,i] = v1
        
        S=np.zeros([m,n])
    
        for i in range(r):
            S[i,i]=s[i]
            
        
            
        return v2,S,u.conjugate().T
    






if __name__ == '__main__':
    
    m,n=2,3
  
    a=produce(m, n)
    
    print("produce a ",m,"*",n," matrix randomly")
    print(a)
    print("\n")

    U,S,V=SVD(a)
    print("this matrix can be write as the product of three matrixs as follow")
    print(U)
    print("\n")
    print(S)
    print("\n")
    print(V)
    print("\n")
    
    print("the product of three matrixs is the same as a")
    print(np.dot(U,S).dot(V))


//print matrix aa
void print(int n,int m,complex<double>* aa)
{
    for(int j = 0;j<n;j++)
    {
        for(int i = 0;i<m;i++)
        {
            cout << aa[j*m+i] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}
