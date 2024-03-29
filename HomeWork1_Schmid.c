#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define n 5 //n is the dimension of the space
#define eps 0.01  //eps is the permission error
typedef double MATRIX[20][20];

void Schmidt_orthogonalization(MATRIX vec_1,int dimension);
void Matrix_print(MATRIX vec,int dimension);
void Matrix_random_generation(MATRIX vec,int dimension);
void Matrix_tran(MATRIX vec_1,MATRIX vec_2,int dimension);
void Matrix_product(MATRIX vec_1,MATRIX vec_2,MATRIX vec_result,int dimension);
double Matrix_scalar_product(MATRIX vector_group_1,MATRIX vector_group_2,int dimension,int a,int b);
double Matrix_det(MATRIX vec, int dimension);
int exam(MATRIX vec,int dimension);


int main()
{

    printf("\n the dimension of the matrix is n=%d,and n should be a natural number less than 20",n);
    MATRIX vec={0};
    Matrix_random_generation(vec,n);// generate a random matrix with n dimensions
    printf("\nthe initial matrix is: \n");
    Matrix_print(vec,n);

    //exam weather the matrix can be expanded into space Rn
    double a=Matrix_det(vec,n);
    if(a==0)  printf("\nThis matrix cannot be expanded into space Rn, initialization failed\n");
    else      printf("\n by examing we can ensure that the matrix can be expanded into space Rn. \n");


    Schmidt_orthogonalization(vec,n);// make vec an orthonormal matrix by the schmidt's way
    printf("\nthe orthonormal matrix is: \n");
    Matrix_print(vec,n);


    MATRIX vec_tran={0};//exam the result
    MATRIX vec_exam={0};
    Matrix_tran(vec_tran,vec,n);
    Matrix_product(vec,vec_tran,vec_exam,n); 
    if(exam(vec_exam,n)==1)  printf("\n by examing we can ensure the matrix is orthinirmal \n"); 
    if(exam(vec_exam,n)==0)  printf("\n can not pass the examination! \n");

    double x[20];
    double u[20];
    for(int i=0;i<=n;i++)
    {
        x[i]=drand48();
    }

    printf("\nwe can generate a vectory by random function\n");
    for(int i=0;i<n;i++)
    {
        printf("%f \n",x[i]);
    }
    //consider that matrix_tran is the inverse matrix of initial matrix,we can obtain the new vectory by matrix_tran*x

    double sum=0;
    for(int i=0;i<n;i++)
    {
        sum=0;
        for(int j=0;j<n;j++)
        {
            sum+=x[j]*vec_tran[i][j];
        }
        u[i]=sum;
    }

    printf("\nthe vectory can be expressed as followed based on the matrix\n");
    for(int i=0;i<n;i++)
    {
        printf("%f \n",u[i]);
    }


    return 0;
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

double Matrix_det(MATRIX vec, int dimension)
{
    double b[20][20] = {{0}}; /*定义数组b并初始化*/
    int i = 0, j = 0; /*i,j为行与列,sum为行列式的值*/
    double sum=0;
    int x = 0, c = 0, p = 0; /*用x判断加与减,c,p为中间变量*/
    if (dimension == 1)
        return vec[0][0];
    for (i = 0; i < dimension; i++) /*此处大循环实现将余子式存入数组b中*/
    {
        for (c = 0; c < dimension - 1; c++)
        {
            for (j = 0; j < dimension - 1; j++)
            {
                if (c < i) { /*借助c判断每行的移动方法*/
                    p = 0; /*当p=0时,行列式只向左移,即消去对应的第一列的数*/
                }
                else { /*否则行列式左移后再上移*/
                    p = 1;
                }
                b[c][j] = vec[c + p][j + 1];
            }
        }
        if (i % 2 == 0) { /*i+j（此时j=0,故只考虑i）为偶数,加法预算*/
            x = 1;
        }
        else { /*i+j为奇数,减法运算*/
            x = (-1);
        }
        sum += vec[i][0] * Matrix_det(b,dimension - 1) * x; /*计算行列式的值*/
    }
    return sum; /*将值返回*/
}
