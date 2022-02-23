#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define n 3  //n is the dimension of the space
#define eps 0.01  //eps is the permission error
typedef double MATRIX[20][20];

void Schmidt_orthogonalization(MATRIX vec_1,int dimension);
void Matrix_print(MATRIX vec,int dimension);
void Matrix_random_generation(MATRIX vec,int dimension);
void Matrix_tran(MATRIX vec_1,MATRIX vec_2,int dimension);
void Matrix_product(MATRIX vec_1,MATRIX vec_2,MATRIX vec_result,int dimension);
double Matrix_scalar_product(MATRIX vector_group_1,MATRIX vector_group_2,int dimension,int a,int b);
int exam(MATRIX vec,int dimension);

int main()
{
    MATRIX vec={0};
    Matrix_random_generation(vec,n);// generate a random matrix with n dimensions
    printf("\nthe initial matrix is: \n");
    Matrix_print(vec,n);


    Schmidt_orthogonalization(vec,n);// make vec an orthonormal matrix by the schmidt's way
    printf("\nthe orthonormal matrix is: \n");
    Matrix_print(vec,n);


    MATRIX vec_tran={0};//exam the result
    MATRIX vec_exam={0};
    Matrix_tran(vec_tran,vec,n);
    Matrix_product(vec,vec_tran,vec_exam,n); 
    if(exam(vec_exam,n)==1)  printf("\n by examing we can ensure the matrix is orthinirmal \n"); 
    if(exam(vec_exam,n)==0)  printf("\n can not pass the examination! \n");
}


void Schmidt_orthogonalization(MATRIX vec_1,int dimension)
{
    dimension=n;
    MATRIX vec_2={0};
    double xv,vv,sum;

    for(int i=0;i<n;i++)
    {
        vec_2[i][0]=vec_1[i][0];
    }


    for(int i=1;i<n;i++)
    {
 
        for(int j=0;j<n;j++)
        {
            sum=0; 
            for(int k=0;k<i;k++)
            {
                xv=0;
                vv=0;
                for(int l=0;l<n;l++)
                {
                    xv+=vec_1[l][i]*vec_2[l][k];
                    vv+=vec_2[l][k]*vec_2[l][k];
                }
                sum+=vec_2[j][k]*xv/vv;
            }
            vec_2[j][i]=vec_1[j][i]-sum;
        }
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            vec_1[i][j]=vec_2[i][j]/pow(Matrix_scalar_product(vec_2,vec_2,n,j,j),0.5);
        }
    }
}

void Matrix_print(MATRIX vec,int dimension)

{
    for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            printf("%f   ",vec[i][j]);
        }
        printf("\n");
    }
}

void Matrix_random_generation(MATRIX vec,int dimension)
{
        for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            vec[i][j]=drand48();
        }
    }
}

void Matrix_tran(MATRIX vec_1,MATRIX vec_2,int dimension)
{
    for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            vec_1[i][j]=vec_2[j][i];
        }
    }

}

void Matrix_product(MATRIX vec_1,MATRIX vec_2,MATRIX vec_result,int dimension)
{
    for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            double sum=0.0;
            for(int k=0;k<dimension;k++)
            {
                sum+=vec_1[i][k]*vec_2[k][j];
            }
            vec_result[i][j]=sum;
        }
    }
}

double Matrix_scalar_product(MATRIX vector_group_1,MATRIX vector_group_2,int dimension,int a,int b)
{
    double result = 0;
    for(int i = 0;i < dimension;i ++)
        result += vector_group_1[i][a] * vector_group_2[i][b];
    
    return result;
}

int exam(MATRIX vec,int dimension)
{
    double sum1=0;
    double sum2=0;
    int bull=0;
    for(int i=0;i<dimension;i++)
    {
        sum1+=vec[i][i];
    }

    if ((sum1-dimension)>eps)  bull++;

    for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            sum2+=vec[i][j];
        }
    }

    if((sum2-dimension)>eps)  bull++;

    if(bull==0)  return 1;
    else  return 0;
}