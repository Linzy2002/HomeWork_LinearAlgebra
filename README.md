# 1、复数域矩阵运算

     复数域矩阵运算可通过函数库BLSA以及complex实现，由于在编写施密特正交化过程中自己整理了一些函数，就此抛弃于心不忍，因此，本次作业调用了自己编写的matrix.h，以实现矩阵乘法，求矩阵转置复共轭等操作。
     

# 2、生成随机向量

本次作业需要生成两个等长度向量$xi$，$yita$。向量维度为N，作为参数。在程序开头通过#define N  确定N的值。由于向量为随机生成，因此不能保证向量长度相等，需要求向量的范数，调整向量长度。这个过程通过函数$void\quad produce(int \quad n,complex<double>* xi, complex<double>* yita);$

实现。

![](C:\Users\lzy\Desktop\3.png)

![](C:\Users\lzy\Desktop\6.png)

参数N不同，则生成的向量维度不同，矩阵的大小不同

# 3、GIVENS 实现

givens矩阵是确定沟通两个等长度向量的初等旋转变换。

通过$int \quad Givens(int n, complex<double>* xi, complex<double>* yita, complex<double>* G)$

函数实现从向量$xi$，和向量$yita$得到沟通他们的givens矩阵。

<img title="" src="file:///C:/Users/lzy/Desktop/@FT0A__@F3DVQLEK$)}V4{7.png" alt="" width="455">

如上图的矩阵可以实现第i个、第j个基矢张成的平面中的旋转。此函数首先对矢量$xi$进行操作，通过旋转$xi$除第一个分量外全部转为0。则可得到

$Gt_1\vec{xi}=\vec{xi_{0}}$，

$xi_{0}$仅第一个分量不为0.同理，可得到

$Gt_2\vec{yita}=\vec{yita_{0}}$,

$yita_0$仅第一个分量不为0.由于旋转变化保内积，因此

$$
\vec{yita_0}=\vec{xi_0}
$$

可得到$Gt_1\vec{xi}=Gt2\vec{yita}$

初等旋转矩阵为幺正矩阵，由幺正矩阵性质${Gt2}^{\dagger}{Gt2}=I$

我们可以得到沟通两向量的矩阵$G\vec{xi}=({Gt2}^{\dagger}Gt1)\vec{xi}=\vec{yita}$。

$Gt$矩阵可通过简单的数学运算求得，则通过矩阵乘法可得到givens矩阵。

![](C:\Users\lzy\Desktop\8.png)

这是N=3时对应某一组随即向量的G矩阵计算结果



# 4、HOUSEHOLDER实现

householder矩阵为一个代表了反射变换的幺正矩阵。可以表示为

$$
H=I-2\omega \omega^{\dagger}  
\tag1
$$

其中$\omega$为垂直反射轴的单位向量。可通过两向量相减得到

<img title="" src="file:///C:/Users/lzy/Desktop/M`PAUG~A5V17@]GD9L%0HCJ.png" alt="" width="318">

由于为复空间中的计算，还需进行如下修正：$ omega[i] = eitheta * xi[i] - yita[i];$

其中$eitheta$代表着相位因子。得到$\omega$矩阵后，将其与自身复共轭相乘，带入$(1)$式中，得到householder矩阵。



因此，得到了计算H矩阵的函数

$int\quad Householder(int\quad n, complex<double>* xi, complex<double>* yita, complex<double>* H)$

![](C:\Users\lzy\Desktop\7.png)

这是N=3时，对应某一组随机向量的H矩阵计算结果

# 5、验证幺正矩阵

幺正矩阵满足性质$UU^{\dagger}=I$，因此只需进行矩阵乘法，将$UU^{\dagger}$与单位阵比较，若为单位阵，则矩阵$U$为幺正矩阵。

![](C:\Users\lzy\Desktop\2.png)





# 6、检验计算结果

程序得到沟通两个等长向量$xi$，$yita$的矩阵，因此检验矩阵是否正确仅需将得到的矩阵作用在向量$xi$上，若能得到$yita$，且矩阵形式满足反射矩阵或旋转矩阵，则证明函数得到了正确的矩阵

![](C:\Users\lzy\Desktop\3.png)

![](C:\Users\lzy\Desktop\4.png)

![](C:\Users\lzy\Desktop\5.png)

这里以第一组向量为例，可以看到，G矩阵和H矩阵作用在$xi$上都得到了$yita$，说明我们得到了正确的旋转矩阵和反射矩阵。
