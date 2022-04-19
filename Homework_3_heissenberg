#include<stdio.h>
#include<iostream>
#include<math.h>

#include "globalconfig.h"
#include "cppblas.h"
#include "cppblasex.h"
#include "Complex.h"
#include "cpplapack.h"



#define L 4 //定义海森堡链的长度


void kron(int m,double* M,int n,double* N,double* Ans);
void print(int n,int m,double* aa);
void matrix_plus(int m,double* M,double* N,double* A);



int main()
{
   

    static double H[(int)pow(2,L)*(int)pow(2,L)]={0};  //H为最终的哈密顿矩阵。
    static double T1[(int)pow(2,L)*(int)pow(2,L)]={0};
    static double T2[(int)pow(2,L)*(int)pow(2,L)]={0};
    static double zeros[(int)pow(2,L)*(int)pow(2,L)]={0};
    static double I2[(int)pow(2,L)*(int)pow(2,L)]={0};
    static double Sz[(int)pow(2,L)*(int)pow(2,L)]={0}; 
    static double Hn[(int)pow(2,L)*(int)pow(2,L)]={0};
    static double Hv[(int)pow(2,L)]={0};
    static double Ht[(int)pow(2,L)]={0};


        
    double I[4]={1,0,0,1};
    double S[16]={0.25,0,0,0,0,-0.25,0.5,0,0,0.5,-0.25,0,0,0,0,0.25};
    double s1[4]={0,1,0,0};
    double s2[4]={0,0,1,0};
    double sz[4]={0.5,0,0,-0.5};


    

    
for(int i=0;i<L-1;i++)
{
    for(int j=0;j<(int)pow(2,L-2)*(int)pow(2,L-2);j++)
    {
        I2[j]=0;
    }

    for(int j=0;j<(int)pow(2,i);j++)
    {
        I2[j*(int)pow(2,i)+j]=1;
    }

    kron((int)pow(2,i),I2,4,S,T1);

    for(int j=0;j<(int)pow(2,L-2)*(int)pow(2,L-2);j++)
    {
        I2[j]=0;
    }

    for(int j=0;j<(int)pow(2,L-i-2);j++)
    {
        I2[j*(int)pow(2,L-i-2)+j]=1;
    }

    kron((int)pow(2,i+2),T1,(int)pow(2,L-i-2),I2,T2);


    matrix_plus((int)pow(2,L),H,T2,H);


}

//最后一个矩阵

    for(int j=0;j<(int)pow(2,L-2)*(int)pow(2,L-2);j++)
    {
        I2[j]=0;
    }

    for(int i=0;i<(int)pow(2,L-2);i++)
    {
        I2[i*(int)pow(2,L-2)+i]=1;
    }

 

    kron(2,s1,(int)pow(2,L-2),I2,T1);
    kron((int)pow(2,L-1),T1,2,s2,T2);
     
    for(int i=0;i<(int)pow(2,L)*(int)pow(2,L);i++)
    {
        T2[i]=0.5*T2[i];
    }
    matrix_plus((int)pow(2,L),H,T2,H);
    
    
    
    kron(2,s2,(int)pow(2,L-2),I2,T1);
    kron((int)pow(2,L-1),T1,2,s1,T2);

    for(int i=0;i<(int)pow(2,L)*(int)pow(2,L);i++)
    {
        T2[i]=0.5*T2[i];
    }
    matrix_plus((int)pow(2,L),H,T2,H);

    
    
    kron(2,sz,(int)pow(2,L-2),I2,T1);
    kron((int)pow(2,L-1),T1,2,sz,T2);
    matrix_plus((int)pow(2,L),H,T2,H);


    //生成Sz矩阵

for(int i=0;i<L;i++)
{
    for(int j=0;j<(int)pow(2,L-1)*(int)pow(2,L-1);j++)
    {
        I2[j]=0;
    }

    for(int j=0;j<(int)pow(2,i);j++)
    {
        I2[j*(int)pow(2,i)+j]=1;
    }

    kron((int)pow(2,i),I2,2,sz,T1);

    for(int j=0;j<(int)pow(2,L-1)*(int)pow(2,L-1);j++)
    {
        I2[j]=0;
    }

    for(int j=0;j<(int)pow(2,L-i-1);j++)
    {
        I2[j*(int)pow(2,L-i-1)+j]=1;
    }

    kron((int)pow(2,i+1),T1,(int)pow(2,L-i-1),I2,T2);


    matrix_plus((int)pow(2,L),Sz,T2,Sz);


}

//选出sz中特征值相同的部分

int number1=0;
int number2=0;
int site1=0;
int site2=0;

double hvt=0;

for(int i=0;i<L+1;i++)
{
    number1=0;
    for(int j=0;j<(int)pow(2,L);j++)
    {
        if (abs(Sz[j*(int)pow(2,L)+j]-(i-L/2))<0.001)
        {   
            Ht[number1]=j;
            number1++;
        }
    }

    for(int j=0;j<number1;j++)
    {
        for(int k=0;k<number1;k++)
        {
            site1=Ht[j];
            site2=Ht[k];
            Hn[j*(number1)+k]=H[site1*(int)pow(2,L)+site2];
        }
    }


    //get the eigenvalue of H

    heev(number1,Hn,number1,Ht,'D');

    for(int j=0;j<number1;j++)
    {
        Hv[number2]=Ht[j];
        number2++;
    }
}



for(int i=0;i<number2;i++)
	{
		for(int j=i+1;j<number2;j++)
		{
			if(Hv[j]<Hv[i])
			{
				hvt=Hv[i];
				Hv[i]=Hv[j];
				Hv[j]=hvt;
			}
		}
	}

    
FILE *fp2;
fp2=fopen("H.dat","w");
for(int i=0;i<(int)pow(2,L);i++)
{
    for(int j=0;j<(int)pow(2,L);j++)
    {
	    fprintf(fp2," %f  ",H[i*(int)pow(2,L)+j]);
    }
    fprintf(fp2,"\n");
}
fclose(fp2);


FILE *fp3;
fp3=fopen("Sz.dat","w");
for(int i=0;i<(int)pow(2,L);i++)
{
    for(int j=0;j<(int)pow(2,L);j++)
    {
	    fprintf(fp3," %f  ",Sz[i*(int)pow(2,L)+j]);
    }
    fprintf(fp3,"\n");
}
fclose(fp3);

FILE *fp4;
fp4=fopen("Hv.dat","w");
for(int i=0;i<number2;i++)
{
	fprintf(fp4," %f  ",Hv[i]);
}
fclose(fp4);

    return 0;
}

//定义直积
void kron(int m,double* M,int n,double* N,double* Ans)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<n;k++)
            {
                for(int l=0;l<n;l++)
                {
                    Ans[(i*n+k)*m*n+(j*n+l)]=M[i*m+j]*N[k*n+l];
                }
            }
        }
    }

}

//打印矩阵
void print(int n,int m,double* aa)
{
    for(int j = 0;j<n;j++){
        for(int i = 0;i<m;i++)
        {
           printf("%f  ",aa[j*m+i]);
        }
        printf("\n");
    }
    printf("\n");
}


//矩阵加法
void matrix_plus(int m,double* M,double* N,double* A)
{
    
    static double temp[(int)pow(2,L)*(int)pow(2,L)]={0};
   
   for(int i=0;i<m*m;i++)
    {
        temp[i]=M[i];
    }

    for(int i=0;i<m*m;i++)
    {
        A[i]=temp[i]+N[i];
    }
}
