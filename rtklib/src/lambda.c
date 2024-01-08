/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    // 调序变换
    j=n-2; k=n-2;   

    // 对第 0,1，...，k-1，k 列进行降相关
    while (j>=0) {
        // 从最后一列开始，各列非对角线元素从上往下依次降相关
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];

        // 检验条件，若不满足检验条件则开始进行调序变换
        // 完成调序变换后重新从最后一列开始进行降相关及排序，k记录最后一次进行过调序变换的列序号
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0; //k表示当前层，从最后一层（n-1）开始计算
    zb[k]=zs[k];//即zn
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);    //四舍五入取整；取整后的数与未取整的数作差；step记录z[k]是四舍还是五入
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {      //如果当前累积目标函数计算值小于当前超椭圆半径
            //情况1：若还未计算至第一层，继续计算累积目标函数值
            if (k!=0) { 
                dist[--k]=newdist;  //记录下当前层的累积目标函数值，dist[k]表示了第k,k+1,...,n-1层的目标函数计算和
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];   //计算Zk，即第k个整数模糊度参数的备选组的中心
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);    //四舍五入取整；取整后的数与未取整的数作差；记录是四舍还是五入
            }
            //情况2：若已经计算至第一层，意味着所有层的累积目标函数值计算完毕
            else {
                //nn为当前候选解数，m为我们需要的固定解数，这里为2，表示需要一个最优解及一个次优解
                //s记录候选解的目标函数值，imax记录之前候选解中的最大目标函数值的坐标
                if (nn<m) { //若候选解数还没满
                    if (nn==0||newdist>s[imax]) imax=nn;    //若当前解的目标函数值比之前最大的目标函数值都大，那么更新imax使s[imax]指向当前解中具有的最大目标函数值
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];  //zn存放所有候选解
                    s[nn++]=newdist;    //s记录当前目标函数值newdist，并加加当前候选解数nn
                }
                else {  //若候选解数已满（即当前zn中已经存了2个候选解）
                    if (newdist<s[imax]) {  //若当前解的目标函数值比s中的最大目标函数值 
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];    //用当前解替换zn中具有较大目标函数值的解
                        s[imax]=newdist;    //用当前解的目标函数值替换s中的最大目标函数值
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;    //更新imax保证imax始终指向s中的最大目标函数值
                    }
                    maxdist=s[imax];    //用当前最大的目标函数值更新超椭圆半径
                }
                //在第一层，取下一个有效的整数模糊度参数进行计算（若zb为5.3，则z取值顺序为5,6,4,7，...）
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        //情况3：如果当前累积目标函数计算值大于当前超椭圆半径
        else {
            if (k==n-1) break;  //如果当前层为第n-1层，意味着后续目标函数各项的计算都会超出超椭圆半径，因此终止搜索
            else {  //若当前层不是第n-1层
                k++;    //退后一层，即从第k层退到第k+1层
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]); //计算退后一层后，当前层的下一个有效备选解
            }
        }
    }
    // 对s中的目标函数值及zn中的候选解进行排序（以s中目标函数值为排序标准，进行升序排序）
    // RTKLIB中最终可以得到一个最优解一个次优解，存在zn中，两解对应的目标函数值，存在s中
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  浮点解数量
*          int    m      I  固定解数量
*          double *a     I  浮点参数向量 (n x 1)
*          double *Q     I  浮点参数协方差阵 (n x n)
*          double *F     O  固定解 (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    // 调用 LD() 首先对浮点协方差阵进行 LD 分解
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        // 调用 reduction() 进行 lambda 降相关性
        /* lambda reduction */
        reduction(n,L,D,Z);

        //  z变换，将双差模糊度进行变换
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        // 调用 search() 进行 mlambda search，将结果存储在 E 和 s 中（整数解）
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            // 逆 Z 变换，将在新空间中固定的模糊度逆变换回双差模糊度空间中，存储在 F 中
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}
