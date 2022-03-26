#include<iostream>
#include<cmath>
#include"matrix.h"
#include<complex>
#include<stdlib.h>
#include<time.h>

#define N 3    //空间维数

using namespace std;

void produce(int n,complex<double>* xi, complex<double>* yita);

int Householder(int n, complex<double>* xi, complex<double>* yita, complex<double>* H);

int Givens(int n, complex<double>* xi, complex<double>* yita, complex<double>* G);

int main()
{
    for(int linzy=0;linzy<3;linzy++)
    {
        cout<<"第"<<linzy+1<<"组向量"<<"\n";
    int n = N; //向量的维度为定义的参数N

    complex<double> xi[n] = {0},yita[n] = {0},r[n] = {0};  //初始化向量
    complex<double> U[n*n] = {0},H[n*n] = {0}, G[n*n] = {0},H_daige[n*n]={0},G_daige[n*n]={0},I[n*n]={0};  //初始化矩阵

    produce(n,xi,yita);//随机生成向量

    Householder(n,xi,yita,H);//通过反射变换生成幺正矩阵

    cout << "H:\n";
    print(n,n,H);   //打印生成的矩阵


    matrix_multiplicate(n,n,n,H,matrix_conjugate(n,n, matrix_transpose(n,n,H)),I);
     cout << "HH_daige=\n";
     print(n,n,I);
    cout<<" HH_daige为单位阵，H为幺正矩阵\n"; //验证幺正矩阵

    matrix_multiplicate(n,n,1,H,xi,r);//将H矩阵作用在xi上

    cout << "\nH*xi:";
    print(1,n,r);

    Givens(n,xi,yita,G);  //通过givens变换得到幺正矩阵
    cout << "G:\n";
    print(n,n,G);

    
    matrix_multiplicate(n,n,n,G,matrix_conjugate(n,n, matrix_transpose(n,n,G)),I);
     cout << "GG_daige=\n";
     print(n,n,I);
    cout<<" GG_daige为单位阵，G为幺正矩阵\n";//验证幺正矩阵

    matrix_multiplicate(n,n,1,G,xi,r);
    cout << "\n G*xi:";
    print(1,n,r);
    }

    return 0;
}


//生成向量
void produce(int n,complex<double>* xi, complex<double>* yita)
{
    double lenth_xi = 0,lenth_yita = 0;
    //srand((unsigned int)time(NULL));
    for(int i = 0;i<n;i++){
        xi[i]  = complex <double> (0.7+0.3*(((double)rand())/RAND_MAX),0.7+0.3*(((double)rand())/RAND_MAX));
        yita[i]  = complex <double> (0.7+0.3*(((double)rand())/RAND_MAX),0.7+0.3*(((double)rand())/RAND_MAX));
        lenth_xi += norm(xi[i]);
        lenth_yita += norm(yita[i]);
    }
    lenth_xi = sqrt(lenth_xi);
    lenth_yita = sqrt(lenth_yita);

    //两向量长度相等
    for(int i = 0;i<n;i++){
        yita[i]  *= lenth_xi/lenth_yita;
    }
    lenth_xi = 0;
    lenth_yita = 0;
    for(int i = 0;i<n;i++){
        lenth_xi += norm(xi[i]);
        lenth_yita += norm(yita[i]);
    }
    cout <<"xi:";
    print(1,n,xi);
    cout <<"yita:";
    print(1,n,yita);
}

int Givens(int n, complex<double>* xi, complex<double>* yita, complex<double>* G)
{
    complex<double> *txi = new complex<double>[n];
    copy(n,xi,txi);
    complex<double> nor = xi[0]/yita[0]; 
    int num = 0;
    
    for (int i = 0; i<n; i++){
        // if (xi[i]/yita[i] == nor) num++;
        G[i*n+i] = complex<double>(1);        
    }

    // if (num == n){
    //     for (int i = 0; i<n; i++) G[i*n+i] = nor;    
    // }

    complex<double> *temg = new complex<double>[n*n],*temxi = new complex<double>[n],sin,cos;
   
   
    for (int ii = 1;ii<n;ii++){
        complex<double>* Gt = new complex<double>[n*n];
       
        for (int i = 0; i<n; i++)
        {
            Gt[i*n+i] = complex<double>(1); // 单位阵
        }

        sin = xi[ii]/sqrt(xi[0]*xi[0] + xi[ii]*xi[ii]);
        cos = xi[0]/sqrt(xi[0]*xi[0] + xi[ii]*xi[ii]);   //确定一次旋转的参数

        Gt [0*n+0] = cos;
        Gt [ii*n+ii] = cos;
        Gt [0*n+ii] = sin;
        Gt [ii*n+0] = -sin;

        matrix_multiplicate(n,n,n,Gt,G,temg);       //将一次旋转作用至原矩阵上
        copy(n*n,temg,G);
        matrix_multiplicate(n,n,1,Gt,xi,temxi);
        copy(n,temxi,xi);
    }

    complex<double>* G2 = new complex<double>[n*n];
    
    for (int i = 0; i<n; i++){
        G2[i*n+i] = complex<double>(1); // 单位阵
    }

    for (int ii = 1;ii<n;ii++){
        complex<double>* Gt = new complex<double>[n*n];

        for (int i = 0; i<n; i++){
            Gt[i*n+i] = complex<double>(1); // 单位阵
        }

        sin = yita[ii]/sqrt(yita[0]*yita[0] + yita[ii]*yita[ii]);
        cos = yita[0]/sqrt(yita[0]*yita[0] + yita[ii]*yita[ii]);

        Gt [0*n+0] = cos; 
        Gt [ii*n+ii] = cos;
        Gt [0*n+ii] = sin;
        Gt [ii*n+0] = -sin;

        matrix_multiplicate(n,n,n,Gt,G2,temg);
        copy(n*n,temg,G2);
        matrix_multiplicate(n,n,1,Gt,yita,temxi);
        copy(n,temxi,yita);
        
    }
    
    scalar_multiplicate(yita[0]/xi[0],n,n,G);
    matrix_multiplicate(n,n,n,matrix_transpose(n,n,G2),G,temg);
    copy(n*n,temg,G);
    copy(n,txi,xi);
    return 0;
}

int Householder(int n, complex<double>* xi, complex<double>* yita, complex<double>* H)
{

    complex<double> nor = xi[0]/yita[0];
    
    //int num = 0;
    // for (int i = 0; i<n; i++){          检验线性无关性  没用              
    //     if (xi[i]/yita[i] == nor) num++;   
    // }
    // if (num == n){
    //     for (int i = 0; i<n; i++) H[i*n+i] = nor;
    //     return 0;
    // }

    complex<double> eitheta = (0),tem[1] = {0};
    matrix_multiplicate(1,n,1,matrix_conjugate(n,1,xi),yita,tem);
    eitheta = tem[0]/sqrt(norm(tem[0]));

    complex<double> omega[n] = {0},omega2[n] = {0};
    double lenth = 0;

    for (int i = 0; i<n; i++)      //求omega
    {
        omega[i] = eitheta * xi[i] - yita[i];
        lenth += norm(omega[i]);
    }

    lenth = sqrt(lenth);
    for (int i = 0; i<n; i++)
    {
        omega[i] /= lenth;
        omega2[i] = conj(omega[i]);
    }


    matrix_multiplicate(n,1,n,omega,omega2,H);
    scalar_multiplicate(complex<double>(2,0),n,n,H);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            H[i*n+j] =  complex <double>(-1)*H[i*n+j];
        }
        H[i*n+i] +=  complex <double>(1,0);
    }
    scalar_multiplicate(eitheta,n,n,H);
    return 0; 
}