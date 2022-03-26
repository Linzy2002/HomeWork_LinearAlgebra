#include<iostream>
#include<cmath>
#include<complex>

using namespace std;



void print(int n,int m,complex<double>* aa){
    for(int j = 0;j<n;j++){
        for(int i = 0;i<m;i++){
            cout << aa[j*m+i] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void matrix_multiplicate(int m,int n,int o,complex<double>* aa,complex<double>*bb,complex<double>* cc){
    for (int i = 0; i < m; i++){
        for (int k = 0; k < o; k++){
            cc[i*o+k] = complex <double>(0,0);
            for (int j = 0; j < n; j++){
                cc[i*o+k] += aa[i*n+j] * bb[j*o+k];
            }
        }
    }
}

void scalar_multiplicate(complex<double> a,int m,int n,complex<double>* aa){
    for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                aa[i*m+j] *= a;
        }
    }
}

complex<double>* matrix_subtract(int m,int n,complex<double>* aa,complex<double>*bb){
    complex<double>* cc = new complex<double>[m*n];
    for (int i = 0 ; i<(n*m) ; i++) cc[i] = aa[i] - bb[i];
    return cc;
}

complex<double>* matrix_conjugate(int m,int n,complex<double>* aa){
    complex<double>* cc = new complex<double>[m*n];
    for (int i = 0 ; i<(n*m) ; i++) cc[i] = conj(aa[i]);
    return cc;
}

void copy(int num,complex<double>* aa,complex<double>*bb){
    for (int i = 0; i<num ; i++) bb[i] = aa[i];
}

complex<double>*  matrix_transpose(int m,int n,complex<double>* aa){
    complex<double>* cc = new complex<double>[m*n];
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            cc[j*n + i] = aa[i*m + j];
        }
    }
    return cc;
}