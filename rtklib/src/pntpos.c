/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#if 0 /* enable GPS-QZS time offset estimation */
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

// 判断导航系统伪距测量值的误差 (E6.24) 的第一项
/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    // 两个三目运算符，如果是 GLONASS 误差因子为 1.5，如果是 SBAS 误差因子值为 3，其它 误差因子值为 1
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

    // 如果高度角过小，设为5°
    if (el<MIN_EL) el=MIN_EL;

    // sin(el)没开方，疑似bug，不过影响不大
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));

    // 消电离层组合，方差*3*3
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}

// 获取 TGD 群波延迟参数
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        // 从 NAV 中找对应的星历
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}

// 调用testsnr()，根据接收机高度角和信号频率来检测该信号是否可用
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    if (testsnr(0,0,azel[1],obs->SNR[0]*SNR_UNIT,&opt->snrmask)) {
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        if (testsnr(0,1,azel[1],obs->SNR[1]*SNR_UNIT,&opt->snrmask)) return 0;
    }
    return 1;
}
/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    double P1,P2,gamma,b1,b2;
    int sat,sys;
    
    sat=obs->sat;
    sys=satsys(sat,NULL);
    P1=obs->P[0];
    P2=obs->P[1];
    *var=0.0;
    
    // 如果没有足够伪距观测值，直接返回
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    
    // 从 nav->cbias 取 DCB 数据
    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }

    // 如果是消电离层组合，将 C1、C2 伪距做 DCB 改正，加上 P1_C1、P2_C2 归化到 P1、P2，计算得到消电离层观测值 PC
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
            gamma=SQR(FREQ1/FREQ2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b */
            gamma=SQR(FREQ1/FREQ7);
            if (getseleph(SYS_GAL)) { /* F/NAV */
                P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQ1)/FREQ2_CMP);
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            b2=gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQ5/FREQ9);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }

    // 如果单频，将 C1 伪距做 DCB 改正，归化到 P1
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
            else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQ9/FREQ5);
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}

// 计算电离层改正量的入口函数，调用 ionmodel、sbsioncorr、iontec 实现
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC) {        // 广播星历模型
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {        // SBAS 校正
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {         // TEC 格网模型
        return iontec(time,nav,pos,azel,1,ion,var);
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {   // QZSS广播星历模型
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    
    // 有电离层改正，方差为0，无电离层改正，方差 5*5
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}

// 计算对流层改正量的入口函数，调用 tropmodel、sbstropcorr 来实现
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   信号发射时间
*          nav_t  *nav      I   NAV 导航数据
*          double *pos      I   接收机纬经高 {lat,lon,h} (rad|m)
*          double *azel     I   方位角、高度角 {az,el} (rad)
*          int    tropopt   I   对流层处理选项 (TROPOPT_???)
*          double *trp      O   对流层延迟 (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* Saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}

// 计算当前迭代的伪距残差 v、几何矩阵 H、伪距残差的方差 var、所有观测卫星的方位角和仰角 azel、
// 定位时有效性 vsat、定位后伪距残差 resp、参与定位的卫星个数 ns 和方程个数 nv。
/* pseudorange residuals -----------------------------------------------------
 * args:int      iter      I   迭代次数，在estpos()里迭代调用，第i次迭代就传i
 *      obsd_t   *obs      I   观测量数据
 *      int      n         I   观测量数据的数量
 *      double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *dts      I   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *      double   *vare     I   卫星位置和钟差的协方差 (m^2)
 *      int      *svh      I   卫星健康标志 (-1:correction not available)
 *      nav_t    *nav      I   导航数据
 *      double   *x        I   本次迭代开始之前的定位值,7*1,前3个是本次迭代开始之前的定位值，第4个是钟差，后三个分别是gps系统与glonass、galileo、bds系统的钟差。
 *      prcopt_t *opt      I   处理过程选项
 *      double   *v        O   定位方程的右端部分，伪距残差
 *      double   *H        O   定位方程中的几何矩阵
 *      double   *var      O   参与定位的伪距残差的方差
 *      double   *azel     O   对于当前定位值，所有观测卫星的 {方位角、高度角} (2*n)
 *      int      *vsat     O   所有观测卫星在当前定位时是否有效 (1*n)
 *      double   *resp     O   所有观测卫星的伪距残差，(P-(r+c*dtr-c*dts+I+T)) (1*n)
 *      int      *ns       O   参与定位的卫星的个数
 * return : 定位方程组的方程个数 
 -------------------------------------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns)
{
    gtime_t time;
    double r,freq,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P;
    int i,j,nv=0,sat,sys,mask[NX-3]={0};
    
    trace(3,"resprng : n=%d\n",n);
    
    // 将之前得到的定位解信息赋值给 rr 和 dtr 数组，以进行关于当前解的伪距残差的相关计算
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    
    // 遍历当前历元所有观测值 obs[]  
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;   
        time=obs[i].time;       // time 赋值 OBS 的时间
        sat=obs[i].sat;         // sat 赋值 OBS 的卫星

        // 调用satsys()函数，验证卫星编号是否合理及其所属的导航系统
        if (!(sys=satsys(sat,NULL))) continue;
        
        // 检测当前观测卫星是否和下一个相邻数据重复；重复则不处理这一条，去处理下一条
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(2,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++;
            continue;
        }

        // 处理选项中事先指定定位时排除哪些导航系统或卫星，调用 satexclude() 完成
        /* excluded satellite? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        // 调用 geodist() 函数，计算卫星和当前接收机位置之间的几何距离 r 和接收机到卫星方向的观测矢量
        // 然后检验几何距离是否 >0，此函数中会进行地球自转影响的校正（Sagnac效应）
        /* geometric distance */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        
        if (iter>0) {
            // 调用 satazel() 函数，计算在接收机位置处的站心坐标系中卫星的方位角和仰角；若仰角低于截断值，不处理此数据
            /* test elevation mask */
            if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
            
            // 调用 snrmask()->testsnr()，根据接收机高度角和信号频率来检测该信号是否可用
            /* test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt)) continue;
            
            // 调用 ionocorr() 函数，计算电离层延时 I，所得的电离层延时是建立在 L1 信号上的，
            // 当使用其它频率信号时，依据所用信号频组中第一个频率的波长与 L1 波长的关系，对上一步得到的电离层延时进行修正
            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            dion*=SQR(FREQ1/freq);      // 电离层改正量
            vion*=SQR(FREQ1/freq);      // 电离层改正误差
            
            // 调用 tropcorr() 计算对流层延时 T
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }

        // 调用 tropcorr() 计算对流层延时 T
        /* psendorange with code bias correction */
        if ((P=prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        // 伪距残差 (P-(r+c*dtr-c*dts+I+T))   (E.6.21)
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        // 构建设计矩阵 H 单位向量的反，前 3 行为中计算得到的视线向，第 4 行为 1，其它行为 0
        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }

        // 处理不同系统（GLO、GAL、CMP、IRN、QZS）与 GPS 之间的时间偏差 (ISB)
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        // 调用 varerr() 计算此时的导航系统误差，然后累加计算用户测距误差(URE)
        /* variance of pseudorange error */
        var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }

    // 为了防止不满秩的情况，把矩阵 H 补满秩了
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }

    // 返回值 v 和 resp 的主要区别在于长度不一致， v 是需要参与定位方程组的解算的，维度为 nv*1；
    // 而 resp 仅表示所有观测卫星的伪距残余，维度为 n*1，对于没有参与定位的卫星，该值为 0
    return nv;
}

// 对定位结果进行卡方检验和 GDOP 检验，失败返回 0，成功返回 1
/* validate solution ---------------------------------------------------------
 * args:const double *azel  方位角、高度角
 *      const int *vsat     观测卫星在当前定位时是否有效 (1*n)
 *      int n               
 *      const prcopt_t *opt 处理选项
 *      const double *v     
 *      int nv              观测值数
 *      int nx              待估计参数数
 *      char *msg           错误消息
 -----------------------------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    // 对残差进行卡方检验，观测值数要大于待估计参数数
    /* Chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {    // (E.6.33)  nv-nx-1:多余观测数、chisqr:卡方值表
        sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        return 0;
    }

    // GDOP 检验
    /* large GDOP check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {     // (E.6.34)
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}

// 通过伪距实现绝对定位，计算出接收机的位置和钟差，顺带返回实现定位后每颗卫星的{方位角、仰角}、定位时有效性、定位后伪距残差。
/* estimate receiver position ------------------------------------------------
 * args:obsd_t   *obs      I   观测量数据
 *      int      n         I   观测量数据的数量
 *      double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *vare     I   卫星位置和钟差的协方差 (m^2)
 *      int      *svh      I   卫星健康标志 (-1:correction not available)
 *      nav_t    *nav      I   导航数据
 *      prcopt_t *opt      I   处理过程选项
 *      prcopt_t *opt      I   处理过程选项
 *      sol_t    *sol      IO  结果结构体
 *      double   *azel     IO  方位角和俯仰角 (rad)
 *      int      *vsat     IO  卫星在定位时是否有效
 *      double   *resp     IO  定位后伪距残差 (P-(r+c*dtr-c*dts+I+T))
 *      char     *msg      O   错误消息
 * return:1 表示成功，0 表示出错
 -------------------------------------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
    // 如果是第一次定位，即输入的 sol 为空，则 x 初值为 0
    // 如果之前有过定位，则可以将上一历元的定位值 sol->rr[] 作为该历元定位的初始值
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    // 开始进行单历元内最小二乘迭代定位计算
    // 解算成功的条件是最小二乘解算出参数增量 dx 的模长小于 1E-4，并且 valsol() 判定结果合格
    // 在迭代次数内成功即停止迭代输出结果，超出迭代次数解算还不成功则认为当前历元解算失败
    for (i=0;i<MAXITR;i++) {

        // 首先调用 rescode() 函数，计算当前迭代的伪距残差 v、设计矩阵矩阵 H
        // 伪距残差的方差 var、所有观测卫星的方位角和仰角 azel、定位时有效性 vsat、
        // 定位后伪距残差 resp、参与定位的卫星个数 ns 和方程个数 nv
        /* pseudorange residuals (m) */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);
        
        // 确定方程组中方程的个数要大于未知数的个数
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }

        // 以伪距残差的标准差的倒数作为权重，对 H 和 v 分别左乘权重对角阵，得到加权之后的 H 和 v
        /* weighted by Std */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);   // 这里的权重值是对角阵，也就是认为不同测量值的误差相互独立
            v[j]/=sig;  
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }

        // 调用 lsq() 函数，最小二乘估计得到当前 x 的修改量 dx 和定位误差协方差矩阵中的权系数阵 Q
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }

        // 将最小二乘求得的修正量 dx 反馈到当前 x 值中，得到更新之后的参数 x 值
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }

        // 如果求得的修改量 dx 小于截断因子(目前是1E-4)，则将 x[j] 作为最终的定位结果，
        // 对 sol 的相应参数赋值，之后再调用 valsol() 调用 valsol() 对定位结果进行卡方检验和 GDOP 检验
        // 否则，进行下一次迭代解算
        if (norm(dx,NX)<1E-4) {
            sol->type=0;

            // 解方程时的钟差单位是 m，是乘以了光速之后的，解出结果后赋给 sol->dtr 时再除以光速
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            // 调用 valsol() 对定位结果进行卡方检验和 GDOP 检验
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    return 0;
}

//RAIM：接收机自主完好性检测与故障检测排除，通过每次排除一颗卫星计算，
//      通过选取残差最小的一次作为最终解算结果，此时对应的卫星就是故障卫星
/* RAIM FDE (failure detection and exclution) -------------------------------
 * arge:const obsd_t *obs       OBS 观测数据
 *      int n                   观测数据的数量
 *      const double *rs        卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      const double *dts       卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *      const double *vare      卫星位置和钟差的协方差 (m^2)
 *      const int *svh          卫星健康标志 (-1:correction not available)
 *      const nav_t *nav        导航数据
 *      const prcopt_t *opt     处理过程选项
 *      sol_t *sol              结果结构体
 *      double *azel            方位角和俯仰角 (rad)
 *      int *vsat               卫星在定位时是否有效
 *      double *resp            观测卫星的伪距残差，(P-(r+c*dtr-c*dts+I+T)) (1*n)
 *      char *msg               错误信息
 -----------------------------------------------------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    // 三层 for 循环
    // i 表示最外面的大循环，每次将将第 i 颗卫星舍弃不用，这是通过 if (j==i) continue实现的
    // j 表示剩余使用的卫星的循环，每次进行相应数据的赋值
    // k 表示参与定位的卫星的循环，与 j 一起使用
    for (i=0;i<n;i++) {
        
        // 舍弃第 i 颗卫星后，将剩下卫星的数据复制到一起
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }

        // 调用 estpos() 使用剩下卫星进行定位解算
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }

        // 累加得到当前卫星实现定位后的伪距残差平方和与可用卫星数目
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }

        // 如果 nvsat<5，则说明当前卫星数目过少，无法进行 RAIM_FDE 操作
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }

        // 计算伪距残差平方和的标准差 rms_e
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        // 如果伪距残差平方和的标准差 rms_e>rms，继续下一次循环，排除下一颗卫星
        if (rms_e>rms) continue;
        
        // 如果小于 rms，则说明当前定位结果更合理，将 stat 置为 1
        // 重新更新 sol、azel、vsat(当前被舍弃的卫星，此值置为0)、resp等值，并将当前的 rms_e更新到 rms 中
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }

    // 如果 stat不为 0，则说明在弃用卫星的前提下有更好的解出现，输出信息，指出弃用了哪颗卫星
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}

//计算定速方程组左边的几何矩阵和右端的速度残余，返回定速时所使用的卫星数目
/* range rate residuals ------------------------------------------------------
 * args:obsd_t   *obs      I   观测量数据
 *      int      n         I   观测量数据的数量
 *      double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *dts      I   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *      nav_t    *nav      I   导航数据
 *      double   *rr       I   接收机位置和速度，长度为6，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *x        I   本次迭代开始之前的定速值，长度为4，{vx,vy,vz,drift}
 *      double   *azel     IO  方位角和俯仰角 (rad)
 *      int      *vsat     I   卫星在定速时是否有效
 *      double   *v        O   定速方程的右端部分，速度残差
 *      double   *H        O   定速方程中的几何矩阵
 * return : 定速时所使用的卫星数目
 -------------------------------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H)
{
    double freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    int i,j,nv=0;
    
    trace(3,"resdop  : n=%d\n",n);
    
    // 调用 ecef2pos() 函数，将接收机位置由 ECEF 转换为大地坐标系
    // 调用 xyz2enu() 函数，计算到ENU坐标转换矩阵 E
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        // 去除在定速时不可用的卫星
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }

        // 计算当前接收机位置下 ENU 中的视线向量 e，然后转换得到 ECEF 中视线向量的值
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);
        
        // 计算 ECEF 中卫星相对于接收机的速度，卫星速度rs[j+3+i*6]-传入的定速初值x[j]
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }

        // 计算考虑了地球自转的用户和卫星之间的几何距离变化率
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        // 计算伪距残差（新息向量） v
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        // 构建左端项的设计矩阵 H
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;   // 将观测方程数增 1
    }
    return nv;
}

/* estimate receiver velocity ------------------------------------------------
 * args   : obsd_t *obs      I   OBS观测数据
 *          int      n         I   观测数据的数量
 *          double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *          double   *dts      I   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *          nav_t    *nav      I   导航数据
 *          prcopt_t *opt      I   处理过程选项
 *          sol_t    *sol      IO  结果结构体
 *          double   *azel     IO  方位角和俯仰角 (rad)
 *          int      *vsat     IO  卫星在定位时是否有效
 * return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    // 这里不像定位时，初始值可能为上一历元的位置(从 sol 中读取初始值)，这里定速的初始值直接给定为 0
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv;
    
    trace(3,"estvel  : n=%d\n",n);
    
    v=mat(n,1); H=mat(4,n);
    
    // for 循环进行迭代解算
    for (i=0;i<MAXITR;i++) {
        
        // 调用 resdop() 计算定速方程设计矩阵 H 和新息向量 v，返回定速时所使用的卫星数目
        /* range rate residuals (m/s) */
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) {
            break;
        }

        // 调用最小二乘法 lsq 函数，解出{速度、加速度}的改正量 dx，频漂计算了，但后面并没有存
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        // 将最小二乘计算出的修正量 dx 反馈到 x 中
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        //  检查当前计算出的改正量的绝对值是否小于 1E-6
        //  是：则说明当前解已经很接近真实值了，将接收机三个方向上的速度存入到 sol->rr 中
        //  否：进行下一次迭代解算
        if (norm(dx,4)<1E-6) {
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data                OBS 观测数据
*          int    n         I   number of observation data      OBS 数
*          nav_t  *nav      I   navigation data                 NAV 导航电文数据
*          prcopt_t *opt    I   processing options              处理过程选项
*          sol_t  *sol      IO  solution                        结果
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)     方位角和俯仰角
*          ssat_t *ssat     IO  satellite status              (NULL: no output)     卫星状态
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->stat=SOLQ_NONE;

    // 检验观测值数是否大于0
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }

    // sol->time 赋值第一个观测值的时间
    sol->time=obs[0].time;
    msg[0]='\0';
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    // 如果处理选项不是 SPP，电离层矫正选 Klobuchar 广播星历模型，对流层矫正采用 Saastmoinen 模型
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }

    // 调用 satposs() 计算卫星位置、速度和钟差
    /* satellite positons, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    // 用伪距位置估计，加权最小二乘，其中会调用 valsol() 进行卡方检验和 GDOP 检验
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    // estpos 中 valsol() 检验失败，即位置估计失败，
    // 会调用 raim_fde() 进行 RAIM-FDE 接收机自主完好性监测重新估计，
    // 前提是卫星数 >6、对应参数解算设置 opt->posopt[4]=1
    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }
    
    // 调用 estvel() 用多普勒估计速度，分米级精度，比伪距准的多。
    // 有 BUG，因为读文件的时候多普勒的正负错了
    /* estimate receiver velocity with Doppler */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }

    // 存入方位角和俯仰角
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    // 赋值卫星状态结构体 ssat
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            ssat[i].snr[0]=0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
