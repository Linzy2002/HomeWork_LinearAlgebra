#include<iostream>
#include<cmath>
#include<complex>
#include<stdlib.h>
#include<time.h>

#define M 3    
#define N 3

using namespace std;


void produce(int m,int n,complex<double>* A);
void print(int n,int m,complex<double>* aa);

int main()
{
    complex<double> A[M*N] = {0};
    produce(M,N,A);
    print(M,N,A);

    return 0;
}



//produce matrix A
void produce(int m,int n,complex<double>* A)
{
    //srand((unsigned int)time(NULL));
    for(int i = 0;i<m*n;i++)
    {
       A[i]  = complex <double> (0.2+0.6*(((double)rand())/RAND_MAX),0.2+0.6*(((double)rand())/RAND_MAX));
    }
}


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
