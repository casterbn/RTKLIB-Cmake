/*------------------------------------------------------------------------------
* rtkpos.c : precise positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/12 1.0  new
*           2007/03/13 1.1  add slip detection by LLI flag
*           2007/04/18 1.2  add antenna pcv correction
*                           change rtkpos argin
*           2008/07/18 1.3  refactored
*           2009/01/02 1.4  modify rtk positioning api
*           2009/03/09 1.5  support glonass, gallileo and qzs
*           2009/08/27 1.6  fix bug on numerical exception
*           2009/09/03 1.7  add check of valid satellite number
*                           add check time sync for moving-base
*           2009/11/23 1.8  add api rtkopenstat(),rtkclosestat()
*                           add receiver h/w bias estimation
*                           add solution status output
*           2010/04/04 1.9  support ppp-kinematic and ppp-static modes
*                           support earth tide correction
*                           changed api:
*                               rtkpos()
*           2010/09/07 1.10 add elevation mask to hold ambiguity
*           2012/02/01 1.11 add extended receiver error model
*                           add glonass interchannel bias correction
*                           add slip detectior by L1-L5 gf jump
*                           output snr of rover receiver in residuals
*           2013/03/10 1.12 add otl and pole tides corrections
*           2014/05/26 1.13 support beidou and galileo
*                           add output of gal-gps and bds-gps time offset
*           2014/05/28 1.14 fix bug on memory exception with many sys and freq
*           2014/08/26 1.15 add functino to swap sol-stat file with keywords
*           2014/10/21 1.16 fix bug on beidou amb-res with pos2-bdsarmode=0
*           2014/11/08 1.17 fix bug on ar-degradation by unhealthy satellites
*           2015/03/23 1.18 residuals referenced to reference satellite
*           2015/05/20 1.19 no output solution status file with Q=0
*           2015/07/22 1.20 fix bug on base station position setting
*           2016/07/30 1.21 suppress single solution if !prcopt.outsingle
*                           fix bug on slip detection of backward filter
*           2016/08/20 1.22 fix bug on ddres() function
*           2018/10/10 1.13 support api change of satexclude()
*           2018/12/15 1.14 disable ambiguity resolution for gps-qzss
*           2019/08/19 1.15 fix bug on return value of resamb_LAMBDA()
*           2020/11/30 1.16 support of NavIC/IRNSS in API rtkpos()
*                           add detecting cycle slips by L1-Lx GF phase jump
*                           delete GLONASS IFB correction in ddres()
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0)   /* 初始位置方差 (m^2) */
#define VAR_VEL     SQR(10.0)   /* 初始速度方差 ((m/s)^2) */
#define VAR_ACC     SQR(10.0)   /* 初始加速度方差 ((m/ss)^2) */
#define VAR_HWBIAS  SQR(1.0)    /* 初始 h/w 方差 ((m/MHz)^2) */
#define VAR_GRA     SQR(0.001)  /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15        /* initial zwd (m) */

#define PRN_HWBIAS  1E-6        /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120         /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0        /* max accel for doppler slip detection (m/s^2) */

#define VAR_HOLDAMB 0.001       /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)   /* time sync tolerance for moving-baseline (s) */
                             

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)  // 频率参数数量，一般为解算频率数，消电离层组合为 1
#define NP(opt)     ((opt)->dynamics==0?3:9)                    // 位置参数数量，一般为 3，动力学模式估计速度加速度为 9
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)      // 电离层参数数量，估计 STEC 为卫星数，否则为 0
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))    // 对流层参数数量，不估计对流层为 0，估计对流层为 2 (基准站、流动站各一个)，估计对流层梯度为 6
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)            // GLONASS 频间偏差数量
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))  // 模糊度参数数量
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))           // 除模糊度外参数数量：位置参数NP(opt)+电离层估计参数NI(opt)+对流层参数NT(opt)+GLONASS AR参数NL(opt)
#define NX(opt)     (NR(opt)+NB(opt))                           // 总参数数量

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* 电离层参数下标 (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* 对流层参数下标 (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)      /* 模糊度参数下标 (s:satno,f:freq) */

/* global variables ----------------------------------------------------------*/
static int statlevel=0;          /* rtk status output level (0:off) */
static FILE *fp_stat=NULL;       /* rtk status file pointer */
static char file_stat[1024]="";  /* rtk status file original path */
static gtime_t time_stat={0};    /* rtk status file time */

// 创建解算状态文件，设置输出等级
/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   rtk status file
*          int      level   I   rtk status level (0: off)
* return : status (1:ok,0:error)
* notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
*          The time to replace keywords is based on UTC of CPU time.
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (ns)
*          clk2     : receiver clock bias GLO-GPS (ns)
*          clk3     : receiver clock bias GAL-GPS (ns)
*          clk4     : receiver clock bias BDS-GPS (ns)
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*-----------------------------------------------------------------------------*/
extern int rtkopenstat(const char *file, int level)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];
    
    trace(3,"rtkopenstat: file=%s level=%d\n",file,level);
    
    if (level<=0) return 0;
    
    reppath(file,path,time,"","");
    
    if (!(fp_stat=fopen(path,"w"))) {
        trace(1,"rtkopenstat: file open error path=%s\n",path);
        return 0;
    }
    strcpy(file_stat,file);
    time_stat=time;
    statlevel=level;
    return 1;
}

// 关闭解算状态文件
/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosestat(void)
{
    trace(3,"rtkclosestat:\n");
    
    if (fp_stat) fclose(fp_stat);
    fp_stat=NULL;
    file_stat[0]='\0';
    statlevel=0;
}

// 写解算状态到缓冲区（解算失败不输出）
/* write solution status to buffer -------------------------------------------*/
extern int rtkoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],vela[3]={0},acca[3]={0},xa[3];
    int i,j,week,est,nfreq,nf=NF(&rtk->opt);
    char id[32],*p=buff;
    
    if (rtk->sol.stat<=SOLQ_NONE) {
        return 0;
    }
    // 如果结果状态为 SOLQ_NONE，直接 return 0，不输出
    /* write ppp solution status to buffer */
    if (rtk->opt.mode>=PMODE_PPP_KINEMA) {
        return pppoutstat(rtk,buff);
    }
    est=rtk->opt.mode>=PMODE_DGPS;
    nfreq=est?nf:1;
    tow=time2gpst(rtk->sol.time,&week);
    
    /* receiver position */
    if (est) {
        for (i=0;i<3;i++) xa[i]=i<rtk->na?rtk->xa[i]:0.0;
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->x[0],rtk->x[1],rtk->x[2],xa[0],xa[1],
                   xa[2]);
    }
    else {
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
                   0.0,0.0,0.0);
    }
    /* receiver velocity and acceleration */
    if (est&&rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        if (rtk->na>=6) ecef2enu(pos,rtk->xa+3,vela);
        if (rtk->na>=9) ecef2enu(pos,rtk->xa+6,acca);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],
                   acc[2],vela[0],vela[1],vela[2],acca[0],acca[1],acca[2]);
    }
    else {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->sol.rr+3,vel);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,rtk->sol.dtr[1]*1E9,
               rtk->sol.dtr[2]*1E9,rtk->sol.dtr[3]*1E9);
    
    /* ionospheric parameters */
    if (est&&rtk->opt.ionoopt==IONOOPT_EST) {
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (!ssat->vs) continue;
            satno2id(i+1,id);
            j=II(i+1,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                       rtk->x[j],xa[0]);
        }
    }
    /* tropospheric parameters */
    if (est&&(rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)) {
        for (i=0;i<2;i++) {
            j=IT(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est&&rtk->opt.glomodear==2) {
        for (i=0;i<nfreq;i++) {
            j=IL(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    return (int)(p-buff);
}
// 根据时间拆分解算状态文件
/* swap solution status file -------------------------------------------------*/
static void swapsolstat(void)
{
    gtime_t time=utc2gpst(timeget());   // 获取当前系统时间 time
    char path[1024];
    
    // 如果当前系统时间周内秒与 time_stat 差距小于 1 天，直接 return
    if ((int)(time2gpst(time     ,NULL)/INT_SWAP_STAT)==
        (int)(time2gpst(time_stat,NULL)/INT_SWAP_STAT)) {
        return;
    }
    time_stat=time;
    
    if (!reppath(file_stat,path,time,"","")) {
        return;
    }
    if (fp_stat) fclose(fp_stat);
    
    if (!(fp_stat=fopen(path,"w"))) {
        trace(2,"swapsolstat: file open error path=%s\n",path);
        return;
    }
    trace(3,"swapsolstat: path=%s\n",path);
}
/* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t *rtk)
{
    ssat_t *ssat;
    double tow;
    char buff[MAXSOLMSG+1],id[32];
    int i,j,n,week,nfreq,nf=NF(&rtk->opt);
    
    if (statlevel<=0||!fp_stat||!rtk->sol.stat) return;
    
    trace(3,"outsolstat:\n");
    
    /* swap solution status file */
    swapsolstat();
    
    /* write solution status */
    n=rtkoutstat(rtk,buff); buff[n]='\0';
    
    fputs(buff,fp_stat);
    
    if (rtk->sol.stat==SOLQ_NONE||statlevel<=1) return;
    
    tow=time2gpst(rtk->sol.time,&week);
    nfreq=rtk->opt.mode>=PMODE_DGPS?nf:1;
    
    /* write residuals and status */
    for (i=0;i<MAXSAT;i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        satno2id(i+1,id);
        for (j=0;j<nfreq;j++) {
            fprintf(fp_stat,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.1f,%d,%d,%d,%d,%d,%d\n",
                    week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                    ssat->resp[j],ssat->resc[j],ssat->vsat[j],
                    ssat->snr[j]*SNR_UNIT,ssat->fix[j],ssat->slip[j]&3,
                    ssat->lock[j],ssat->outc[j],ssat->slipc[j],ssat->rejc[j]);
        }
    }
}
/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}

// 计算单差伪距观测值，i 是基准站下标、j 是流动站下标、k 是选择的频率
/* single-differenced observable ---------------------------------------------*/
static double sdobs(const obsd_t *obs, int i, int j, int k)
{
    double pi=(k<NFREQ)?obs[i].L[k]:obs[i].P[k-NFREQ];
    double pj=(k<NFREQ)?obs[j].L[k]:obs[j].P[k-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
// 计算单差频率 GF 组合观测值，i 是基准站下标、j 是流动站下标、k 是选择的频率
/* single-differenced geometry-free linear combination of phase --------------*/
static double gfobs(const obsd_t *obs, int i, int j, int k, const nav_t *nav)
{
    double freq1,freq2,L1,L2;
    
    freq1=sat2freq(obs[i].sat,obs[i].code[0],nav);
    freq2=sat2freq(obs[i].sat,obs[i].code[k],nav);
    L1=sdobs(obs,i,j,0);
    L2=sdobs(obs,i,j,k);
    if (freq1==0.0||freq2==0.0||L1==0.0||L2==0.0) return 0.0;
    return L1*CLIGHT/freq1-L2*CLIGHT/freq2;
}
// 单差量测方差
/* single-differenced measurement error variance -----------------------------*/
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
                     const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double sinel=sin(el);
    int nf=NF(opt);
    
    if (f>=nf) fact=opt->eratio[f-nf];
    if (fact<=0.0) fact=opt->eratio[0];
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    a=fact*opt->err[1];
    b=fact*opt->err[2];
    
    return 2.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;
}
// 求基线长度，流动站坐标减基准站坐标再取模
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}

// 赋值 xi、var、给rtk->x[i]，rtk->P[i,i]，rtk->P[i] 上非对角线元素赋 0
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;               // 赋值 rtk->x[i]
    for (j=0;j<rtk->nx;j++) {   // 遍历 rtk->P[i]，对角线上元素赋方差，对角线外原始赋 0
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

// 选择基准站、流动站间的共视卫星，返回共同观测的卫星个数 ns，输出卫星号列表 sat
// 在接收机观测值中的下标值列表 iu 和在基站观测值中的下标值列表 ir
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
                  const prcopt_t *opt, int *sat, int *iu, int *ir)
{
    int i,j,k=0;
    
    trace(3,"selsat  : nu=%d nr=%d\n",nu,nr);
    
    for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (azel[1+j*2]>=opt->elmin) { /* elevation at base station */
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
            trace(4,"(%2d) sat=%3d iu=%2d ir=%2d\n",k-1,obs[i].sat,i,j);
        }
    }
    return k;
}

// 位置参数时间更新
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos(rtk_t *rtk, double tt)
{
    double *F,*P,*FP,*x,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i,j,*ix,nx;
    
    trace(3,"udpos   : tt=%.3f\n",tt);
    
    // 若为 PMODE_FIXED 模式，直接从选项中取得位置值给 rtk->x，然后返回
    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }

    // 若为第一个历元，用 rtk->sol 中的单点定位解初始化 rtk->x，方差设为 VAR_POS(30*30)
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    // 若为 PMODE_STATIC 模式，不改变协方差阵直接 return
    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;
    
    // 若非 dynamics 模式，用 rtk->sol 中的位置值初始化 rtk->x，方差设为 VAR_POS(30*30)
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        return;
    }

    // 检查位置协方差，若大于阈值 VAR_POS 则用 rtk->sol 中的单点解重置 rtk->x
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx];
    var/=3.0;
    
    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }

    // 生成有效状态量(非0)下标 ix[]
    /* generate valid state index */
    ix=imat(rtk->nx,1);
    for (i=nx=0;i<rtk->nx;i++) {
        if (rtk->x[i]!=0.0&&rtk->P[i+i*rtk->nx]>0.0) ix[nx++]=i;
    }
    if (nx<9) {
        free(ix);
        return;
    }

    // 构建状态转移矩阵 F
    /* state transition of position/velocity/acceleration */
    F=eye(nx); P=mat(nx,nx); FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);
    
    for (i=0;i<6;i++) {
        F[i+(i+3)*nx]=tt;
    }
    for (i=0;i<3;i++) {
        F[i+(i+6)*nx]=SQR(tt)/2.0;
    }
    // rtk->x 赋值参数阵 x[]、rtk->P 赋值方差阵 P[]
    for (i=0;i<nx;i++) {
        x[i]=rtk->x[ix[i]];
        for (j=0;j<nx;j++) {
            P[i+j*nx]=rtk->P[ix[i]+ix[j]*rtk->nx];
        }
    }

    // 用状态转移矩阵对参数和协方差郑进行卡尔曼滤波时间更新
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",nx,1,nx,1.0,F,x,0.0,xp);
    matmul("NN",nx,nx,nx,1.0,F,P,0.0,FP);
    matmul("NT",nx,nx,nx,1.0,FP,F,0.0,P);
    
    for (i=0;i<nx;i++) {
        rtk->x[ix[i]]=xp[i];
        for (j=0;j<nx;j++) {
            rtk->P[ix[i]+ix[j]*rtk->nx]=P[i+j*nx];
        }
    }
    // 给加速度加过程噪声，先存到 Q，再赋值给 rtk->P
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(tt);
    Q[8]=SQR(rtk->opt.prn[4])*fabs(tt);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
    }
    free(ix); free(F); free(P); free(FP); free(x); free(xp);
}

// 电离层参数时间更新
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
    double el,fact;
    int i,j;
    
    trace(3,"udion   : tt=%.3f bl=%.0f ns=%d\n",tt,bl,ns);
    
    // 遍历每颗卫星，如果两个频率，载波中断计数都大于 GAP_RESION(120)，电离层参数设为 0
    for (i=1;i<=MAXSAT;i++) {
        j=II(i,&rtk->opt);
        if (rtk->x[j]!=0.0&&
            rtk->ssat[i-1].outc[0]>GAP_RESION&&rtk->ssat[i-1].outc[1]>GAP_RESION)
            rtk->x[j]=0.0;
    }
    // 遍历共视卫星
    for (i=0;i<ns;i++) {
        j=II(sat[i],&rtk->opt);
        // 如果电离层状态量为 0，状态设为 1E-6，协方差设为 SQR(rtk->opt.std[1]*bl/1E4)，bl为基线长度
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j);
        }
        // 电离层状态量不为 0，加过程噪声 SQR(rtk->opt.prn[1]*bl/1E4*fact)*fabs(tt)
        else {
            /* elevation dependent factor of process noise */
            el=rtk->ssat[sat[i]-1].azel[1];
            fact=cos(el);
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]*bl/1E4*fact)*fabs(tt);
        }
    }
}

// 对流层参数时间更新
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt, double bl)
{
    int i,j,k;
    
    trace(3,"udtrop  : tt=%.3f\n",tt);
    
    // 处理基站、移动站对流层延迟
    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);
        
        // 如果对流层状态量为 0，用 INIT_ZWD、SQR(rtk->opt.std[2])初始化天顶方向对流层延迟
        if (rtk->x[j]==0.0) {
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */
            
            // TROPOPT_ESTG 模式，用 1E-6、VAR_GRA 初始化东向、北向的对流层梯度系数
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) initx(rtk,1E-6,VAR_GRA,++j);
            }
        }
        // 对流层状态量不为 0，增加过程噪声 SQR(rtk->opt.prn[2])*fabs(tt)
        else {
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2])*fabs(tt);
            // TROPOPT_ESTG 模式，给东向、北向的对流层梯度系数加过程噪声`SQR(rtk->opt.prn[2]*0.3)*fabs(tt)`
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    rtk->P[++j*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.3)*fabs(tt);
                }
            }
        }
    }
}
/* temporal update of receiver h/w biases ------------------------------------*/
static void udrcvbias(rtk_t *rtk, double tt)
{
    int i,j;
    
    trace(3,"udrcvbias: tt=%.3f\n",tt);
    
    // 遍历 GLONASS 的频率
    for (i=0;i<NFREQGLO;i++) {
        j=IL(i,&rtk->opt);
        
        // 如果接收机硬件偏移参数为 0，用 1E-6、VAR_HWBIAS 初始化
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,VAR_HWBIAS,j);
        }
        /* hold to fixed solution */
        else if (rtk->nfix>=rtk->opt.minfix&&rtk->sol.ratio>rtk->opt.thresar[0]) {
            initx(rtk,rtk->xa[j],rtk->Pa[j+j*rtk->na],j);
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(PRN_HWBIAS)*fabs(tt);  // 加过程噪声
        }
    }
}

// 根据 LLI 周跳检测
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    uint32_t slip,LLI;
    int f,sat=obs[i].sat;
    
    trace(3,"detslp_ll: i=%d rcv=%d\n",i,rcv);
    
    // nf 是频点个数(1:L1,2:L1+L2,3:L1+L2+L5)
    for (f=0;f<rtk->opt.nf;f++) {
        // bit0 为 1，表明当前历元可能发生了周跳
        if (obs[i].L[f]==0.0||
            fabs(timediff(obs[i].time,rtk->ssat[sat-1].pt[rcv-1][f]))<DTTOL) {
            continue;
        }
        /* restore previous LLI */
        if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
        else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */
        
        /* detect slip by cycle slip flag in LLI */
        if (rtk->tt>=0.0) { /* forward */
            if (obs[i].LLI[f]&1) {
                errmsg(rtk,"slip detected forward  (sat=%2d rcv=%d F=%d LLI=%x)\n",
                       sat,rcv,f+1,obs[i].LLI[f]);
            }
            slip=obs[i].LLI[f];
        }
        else { /* backward */
            // 前一历元 bit0 位为1，表明前一历元发生可能发生了周跳
            if (LLI&1) {
                errmsg(rtk,"slip detected backward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                       sat,rcv,f+1,LLI);
            }
            slip=LLI;
        }

        // (LLI&2)&&!(obs[i].LLI[f]&2)：前一历元 bit1 位为 1 并且当前历元 bit1 位不为 1
        // (!(LLI&2)&&(obs[i].LLI[f]&2)：前一历元 bit1 位为 0 且当前历元 bit1 位为 1
        // 即前后历元的 bit1 位不相同，表明可能发生了半周跳。
        /* detect slip by parity unknown flag transition in LLI */
        if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
            errmsg(rtk,"slip detected half-cyc (sat=%2d rcv=%d F=%d LLI=%x->%x)\n",
                   sat,rcv,f+1,LLI,obs[i].LLI[f]);
            slip|=1;
        }
        /* save current LLI */
        if (rcv==1) setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);
        else        setbitu(&rtk->ssat[sat-1].slip[f],2,2,obs[i].LLI[f]);
        
        /* save slip and half-cycle valid flag */
        rtk->ssat[sat-1].slip[f]|=(uint8_t)slip;
        rtk->ssat[sat-1].half[f]=(obs[i].LLI[f]&2)?0:1;
    }
}

// 利用几何无关组合进行周跳检测
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int i, int j,
                      const nav_t *nav)
{
    int k,sat=obs[i].sat;
    double g0,g1;
    
    trace(3,"detslp_gf: i=%d j=%d\n",i,j);
    
    for (k=1;k<rtk->opt.nf;k++) {
        if ((g1=gfobs(obs,i,j,k,nav))==0.0) return;
         
        g0=rtk->ssat[sat-1].gf[k-1];
        rtk->ssat[sat-1].gf[k-1]=g1;
         
        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            rtk->ssat[sat-1].slip[0]|=1;
            rtk->ssat[sat-1].slip[k]|=1;
            errmsg(rtk,"slip detected GF jump (sat=%2d L1-L%d GF=%.3f %.3f)\n",
                   sat,k+1,g0,g1);
        }
    }
}
/* detect cycle slip by doppler and phase difference -------------------------*/
static void detslp_dop(rtk_t *rtk, const obsd_t *obs, int i, int rcv,
                       const nav_t *nav)
{
#if 0 /* detection with doppler disabled because of clock-jump issue (v.2.3.0) */
    int f,sat=obs[i].sat;
    double tt,dph,dpt,lam,thres;
    
    trace(3,"detslp_dop: i=%d rcv=%d\n",i,rcv);
    
    for (f=0;f<rtk->opt.nf;f++) {
        if (obs[i].L[f]==0.0||obs[i].D[f]==0.0||rtk->ph[rcv-1][sat-1][f]==0.0) {
            continue;
        }
        if (fabs(tt=timediff(obs[i].time,rtk->pt[rcv-1][sat-1][f]))<DTTOL) continue;
        if ((lam=nav->lam[sat-1][f])<=0.0) continue;
        
        /* cycle slip threshold (cycle) */
        thres=MAXACC*tt*tt/2.0/lam+rtk->opt.err[4]*fabs(tt)*4.0;
        
        /* phase difference and doppler x time (cycle) */
        dph=obs[i].L[f]-rtk->ph[rcv-1][sat-1][f];
        dpt=-obs[i].D[f]*tt;
        
        if (fabs(dph-dpt)<=thres) continue;
        
        rtk->slip[sat-1][f]|=1;
        
        errmsg(rtk,"slip detected (sat=%2d rcv=%d L%d=%.3f %.3f thres=%.3f)\n",
               sat,rcv,f+1,dph,dpt,thres);
    }
#endif
}

// 更新载波相位偏移状态（单差整周模糊度）值到 rtk->x 以及更新其误差协方差到 rtk->P
/* temporal update of phase biases -------------------------------------------
 * args:rtk_t    *rtk      IO  rtk控制结构体
 *      double   tt        I   本次更新与上次更新的时间差
 *      obsd_t   *obs      I   观测数据
 *      int      sat       I   接收机和基站共同观测的卫星号列表
 *      int      *iu       I   接收机和基站共同观测的卫星在接收机观测值中的index值列表
 *      int      *ir       I   接收机和基站共同观测的卫星在基站观测值中的index值列表
 *      int      ns        I   接收机和基站共同观测的卫星个数
 *      nav_t    *nav      I   导航数据
 ------------------------------------------------------------------------------*/
static void udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                   const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,cp1,cp2,pr1,pr2,*bias,offset,freqi,freq1,freq2,C1,C2;
    int i,j,k,slip,reset,nf=NF(&rtk->opt);
    
    trace(3,"udbias  : tt=%.3f ns=%d\n",tt,ns);
    
    // 首先是对所有共视星进行周跳检测
    // 遍历每一颗卫星
    for (i=0;i<ns;i++) {
        // 调用 detslp_ll() 根据 LLI 检查接收机和基站观测数据是否有周跳
        /* detect cycle slip by LLI */
        for (k=0;k<rtk->opt.nf;k++) rtk->ssat[sat[i]-1].slip[k]&=0xFC;
        detslp_ll(rtk,obs,iu[i],1);
        detslp_ll(rtk,obs,ir[i],2);
        
        // 调用 detslp_gf() 利用几何无关组合进行周跳检测
        /* detect cycle slip by geometry-free phase jump */
        detslp_gf(rtk,obs,iu[i],ir[i],nav);
        
        // 调用 detslp_dop() 利用多普勒进行周跳检测（但该函数由于时间跳变的原因，暂未使用）
        /* detect cycle slip by doppler and phase difference */
        detslp_dop(rtk,obs,iu[i],1,nav);
        detslp_dop(rtk,obs,ir[i],2,nav);
        
        /* update half-cycle valid flag */
        for (k=0;k<nf;k++) {
            rtk->ssat[sat[i]-1].half[k]=
                !((obs[iu[i]].LLI[k]&2)||(obs[ir[i]].LLI[k]&2));
        }
    }

    // 遍历每一个频率
    for (k=0;k<nf;k++) {

        // 若为 instantaneous 模式，或者卫星载波相位的中断次数rtk->ssat[i-1].outc
        // 大于配置中所设置的最大次数rtk->opt.maxout，重置载波偏移值
        /* reset phase-bias if instantaneous AR or expire obs outage counter */
        for (i=1;i<=MAXSAT;i++) {
            
            reset=++rtk->ssat[i-1].outc[k]>(uint32_t)rtk->opt.maxout;
            
            if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,k,&rtk->opt)]!=0.0) {
                initx(rtk,0.0,0.0,IB(i,k,&rtk->opt));
            }
            else if (reset&&rtk->x[IB(i,k,&rtk->opt)]!=0.0) {
                initx(rtk,0.0,0.0,IB(i,k,&rtk->opt));
                trace(3,"udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
                      i,k+1,rtk->ssat[i-1].outc[k]);
                rtk->ssat[i-1].outc[k]=0;
            }
            if (rtk->opt.modear!=ARMODE_INST&&reset) {
                rtk->ssat[i-1].lock[k]=-rtk->opt.minlock;
            }
        }
        // 对共视星进行循环，对单差相位偏移状态量的协方差阵加入过程噪声 rtk->opt.prn[0]*rtk->opt.prn[0]*fabs(tt)
        // 如果发现有周跳或者异常值，则单差相位偏移状态量重置为 0
        /* reset phase-bias if detecting cycle slip */
        for (i=0;i<ns;i++) {
            j=IB(sat[i],k,&rtk->opt);
            rtk->P[j+j*rtk->nx]+=rtk->opt.prn[0]*rtk->opt.prn[0]*fabs(tt);
            slip=rtk->ssat[sat[i]-1].slip[k];
            if (rtk->opt.ionoopt==IONOOPT_IFLC) slip|=rtk->ssat[sat[i]-1].slip[1];
            if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
            rtk->x[j]=0.0;
            rtk->ssat[sat[i]-1].lock[k]=-rtk->opt.minlock;
        }
        bias=zeros(ns,1);
        
        /* estimate approximate phase-bias by phase - code 
        对共视星进行循环，利用“单差伪距”和“单差载波相位”计算一个“单差相位偏移”估计值，
        来对单差相位偏移状态量进行更新。由于如果忽略伪距误差，那么伪距减去载波相位，
        则应该是载波相位（m）的相位偏移(m)，所以这里计算的是一个大致的相位偏移bias[i]。
        如果配置为无电离层组合，计算则按照无电离层组合的方式来计算。
        最后，仅计算所有有效星（单差相位偏移状态不为0）的offset = sum of (bias - phase-bias)。
        其实就是计算每颗有效星bias[i]与单单差相位偏移状态量的偏差，然后把所有有效星的这个偏差值加起来，
        之后会除以有效星的个数，最终就是求一个偏差平均值rtk->com_bias。
        */
       //对共视星进行循环
        for (i=j=0,offset=0.0;i<ns;i++) {
            
            if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
                cp=sdobs(obs,iu[i],ir[i],k); /* cycle */        // 载波相位单差观测值
                pr=sdobs(obs,iu[i],ir[i],k+NFREQ);              // 伪距单差观测值
                freqi=sat2freq(sat[i],obs[iu[i]].code[k],nav);  // 载波相位频率
                if (cp==0.0||pr==0.0||freqi==0.0) continue;
                
                bias[i]=cp-pr*freqi/CLIGHT;     // 大致的相位偏移
            }
            // 消电离层模式
            else {
                cp1=sdobs(obs,iu[i],ir[i],0);
                cp2=sdobs(obs,iu[i],ir[i],1);
                pr1=sdobs(obs,iu[i],ir[i],NFREQ);
                pr2=sdobs(obs,iu[i],ir[i],NFREQ+1);
                freq1=sat2freq(sat[i],obs[iu[i]].code[0],nav);
                freq2=sat2freq(sat[i],obs[iu[i]].code[1],nav);
                if (cp1==0.0||cp2==0.0||pr1==0.0||pr2==0.0||freq1==0.0||freq2<=0.0) continue;
                
                C1= SQR(freq1)/(SQR(freq1)-SQR(freq2));
                C2=-SQR(freq2)/(SQR(freq1)-SQR(freq2));
                bias[i]=(C1*cp1*CLIGHT/freq1+C2*cp2*CLIGHT/freq2)-(C1*pr1+C2*pr2);
            }
            // 与原来状态量中单差模糊度的偏差的总和
            if (rtk->x[IB(sat[i],k,&rtk->opt)]!=0.0) {
                offset+=bias[i]-rtk->x[IB(sat[i],k,&rtk->opt)]; 
                j++;
            }
        }
        /* correct phase-bias offset to enssure phase-code coherency */
        if (j>0) {
            for (i=1;i<=MAXSAT;i++) {
                if (rtk->x[IB(i,k,&rtk->opt)]!=0.0) rtk->x[IB(i,k,&rtk->opt)]+=offset/j;
            }
        }
        // 利用求得的 bias[i] 来对没有进行初始化的卫星，进行单差相位偏移状态量的初始化
        /* set initial states of phase-bias */
        for (i=0;i<ns;i++) {
            if (bias[i]==0.0||rtk->x[IB(sat[i],k,&rtk->opt)]!=0.0) continue;
            initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],k,&rtk->opt));
        }
        free(bias);
    }
}

// 这个函数主要是对各状态参数进行预测，更新状态值 rtk->x 及其误差协方差 rtk->P
// 状态参数可能包括坐标、速度、加速度、电离层延迟、对流层天顶湿延迟、模糊度参数、GLONASS卫星的接收机硬件延迟，
// 分别对应udpos函数、udion函数、udtrop函数、udbias函数、udrcvbias函数。
/* temporal update of states --------------------------------------------------
 * args:rtk_t    *rtk      IO  rtk控制结构体
 *      obsd_t   *obs      I   观测数据
 *      int      sat       I   接收机和基站共同观测的卫星号列表
 *      int      *iu       I   接收机和基站共同观测的卫星在接收机观测值中的index值列表
 *      int      *ir       I   接收机和基站共同观测的卫星在基站观测值中的index值列表
 *      int      ns        I   接收机和基站共同观测的卫星个数
 *      nav_t    *nav      I   导航数据
 -----------------------------------------------------------------------------*/
static void udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                    const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double tt=rtk->tt,bl,dr[3];
    
    trace(3,"udstate : ns=%d\n",ns);
    
    // 调用 udpos() 根据不同模式更新rtk中的位置、速度、加速度值和协方差
    /* temporal update of position/velocity/acceleration */
    udpos(rtk,tt);
    
    // 电离层模式>=IONOOPT_EST，调用 udion 更新状态 rtk->x 中的电离层参数（MAXSAT个）及其协方差
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        udion(rtk,tt,bl,sat,ns);
    }

    // 若对流层模式>=TROPOPT_EST，调用 udtrop 更新状态 rtk->x 中的对流层参数（2或6个）及其协方差
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        udtrop(rtk,tt,bl);      // 对流层这里也需要用到上面电离层处理时的基线长度 bl
    }
    // 若为 GLONASS AR 模式，调用 udrcvbias() 更新接收机硬件偏移
    /* temporal update of eceiver h/w bias */
    if (rtk->opt.glomodear==2&&(rtk->opt.navsys&SYS_GLO)) {
        udrcvbias(rtk,tt);
    }
    // 若模式>PMODE_DGPS，调用 udbias() 更新载波相位偏移状态值以及其误差协方差
    /* temporal update of phase-bias */
    if (rtk->opt.mode>PMODE_DGPS) {
        udbias(rtk,tt,obs,sat,iu,ir,ns,nav);
    }
}

//计算接收机或基站对某一颗卫星的非差残差（Zero-Difference Residuals）。y = 观测值 - r - dant
/* UD (undifferenced) phase/code residual for satellite ----------------------
 * args:int      base      I   0表示接收机，1表示基站
 *      double   r         I   经过钟差和对流层校正后的几何距离
 *      obsd_t   *obs      I   OBS观测数据
 *      nav_t    *nav      I   NAV导航数据
 *      double   *azel     I   方位角和俯仰角 (rad)
 *      double   *dant     I   接收机天线校正值
 *      prcopt_t *opt      I   处理过程选项
 *      double   *y        O   非差残差
 *      double *freq
 -----------------------------------------------------------------------------*/
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                      const double *azel, const double *dant,
                      const prcopt_t *opt, double *y, double *freq)
{
    double freq1,freq2,C1,C2,dant_if;
    int i,nf=NF(opt);
    
    // 电离层校正模式为 IONOOPT_IFLC 的情况
    if (opt->ionoopt==IONOOPT_IFLC) { /* iono-free linear combination */
        // 获取观测值 obs->code[0]、obs->code[1] 的频率 freq1、freq2
        freq1=sat2freq(obs->sat,obs->code[0],nav);
        freq2=sat2freq(obs->sat,obs->code[1],nav);
        if (freq1==0.0||freq2==0.0) return;
        
        // 调用 testsnr() 检测信噪比
        if (testsnr(base,0,azel[1],obs->SNR[0]*SNR_UNIT,&opt->snrmask)||
            testsnr(base,1,azel[1],obs->SNR[1]*SNR_UNIT,&opt->snrmask)) return;
        
        // 计算消电离层组合系数c1、c2    (E.5.23)
        C1= SQR(freq1)/(SQR(freq1)-SQR(freq2));
        C2=-SQR(freq2)/(SQR(freq1)-SQR(freq2));

        // 计算天线校正值 dant_if
        dant_if=C1*dant[0]+C2*dant[1];
        
        // 计算残差
        if (obs->L[0]!=0.0&&obs->L[1]!=0.0) {   // 载波相位残差  (E.5.22) 
            y[0]=C1*obs->L[0]*CLIGHT/freq1+C2*obs->L[1]*CLIGHT/freq2-r-dant_if;
        }
        if (obs->P[0]!=0.0&&obs->P[1]!=0.0) {   // 伪距残差  (E.5.21) 
            y[1]=C1*obs->P[0]+C2*obs->P[1]-r-dant_if;
        }
        freq[0]=1.0;
    }
    // 电离层校正模式不为 IONOOPT_IFLC 的情况
    else {
        for (i=0;i<nf;i++) {
            if ((freq[i]=sat2freq(obs->sat,obs->code[i],nav))==0.0) continue;            
            /* check SNR mask */    // 检测信噪比
            if (testsnr(base,i,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) {
                continue;
            }
            /* residuals = observable - pseudorange */      // 计算残差
            if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*CLIGHT/freq[i]-r-dant[i];
            if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]               -r-dant[i];
        }
    }
}

// 这个函数主要是计算流动站（base=0）和基准站（base=1）的非差残差，即没有差分的相位/码残差（Zero-Difference Residuals）
// 即卫地距观测值（观测文件内的载波、伪距值）减去卫地距计算值
/* UD (undifferenced) phase/code residuals -----------------------------------
 * args:int      base      I   0表示接收机，1表示基站
 *      obsd_t   *obs      I   OBS观测数据
 *      int      n         I   OBS的数量
 *      double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *dts      I   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *      int      *svh      I   卫星健康标志 (-1:correction not available)
 *      nav_t    *nav      I   NAV导航数据
 *      double   *rr       I   接收机/基站的位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      prcopt_t *opt      I   处理过程选项
 *      int      index     I   0表示接收机，1表示基站，与参数 base 重复了
 *      double   *y        O   相位/码残差
 *      double   *e        O   观测矢量 (ecef)
 *      double   *azel     O   方位角和俯仰角 (rad)
 *      double   *freq     O   载波频率
 * return : (1:ok,0:error)
 -----------------------------------------------------------------------------*/
static int zdres(int base, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const double *var, const int *svh,
                 const nav_t *nav, const double *rr, const prcopt_t *opt,
                 int index, double *y, double *e, double *azel, double *freq)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R};
    int i,nf=NF(opt);
    
    trace(3,"zdres   : n=%d\n",n);
    
    for (i=0;i<n*nf*2;i++) y[i]=0.0;

    // 若没有接收机位置，return 0
    if (norm(rr,3)<=0.0) return 0; /* no receiver position */

    // 接收机位置传给 rr_
    for (i=0;i<3;i++) rr_[i]=rr[i];  
    
    // 若需要地球潮校正，调用 tidedisp() 对 rr_ 进行校正，地球潮包含固体潮、极潮和海潮负荷
    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
                 opt->odisp[base],disp);
        for (i=0;i<3;i++) rr_[i]+=disp[i];
    }
    ecef2pos(rr_,pos);
    
    // 遍历观测量
    for (i=0;i<n;i++) {

        // 根据卫星位置和接收机位置计算卫地距，并进行地球自转改正
        /* compute geometric-range and azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
        if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;
        
        // 排除需要排除的卫星
        /* excluded satellite? */
        if (satexclude(obs[i].sat,var[i],svh[i],opt)) continue;
        
        // 对卫地距进行卫星钟差距离改正r
        /* satellite clock-bias */
        r+=-CLIGHT*dts[i*2];
        
        // 对卫地距进行对流层延迟改正。对流层延迟可分为大概 90% 的干延迟和 10% 的湿延迟，
        // RTKLIB 内利用 Saastamoinen 模型只改正了对流层的干延迟，湿延迟通过随机游走估计的方式进行改正。
        // 此处湿度传0，湿延迟计算出的结果就是 0
        // 当然，建议对这一部分进行算法改进，利用 GPT2_w 等经验模型改正对流层的总延迟
        /* troposphere delay model (hydrostatic) */
        zhd=tropmodel(obs[0].time,pos,zazel,0.0);
        r+=tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;
        
        // 接收机天线相位中心改正，调用 antmodel() 计算校正值 dant（对每一个频率都有一个值）
        /* receiver antenna phase center correction */
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],
                 dant);
        
        // 观测值减经过上述各项改正后的计算值，得到最终的非差残差 v。此处涉及到无电离层组合观测值的求解。
        /* UD phase/code residual for satellite */
        zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*nf*2,freq+i*nf);
    }
    trace(4,"rr_=%.3f %.3f %.3f\n",rr_[0],rr_[1],rr_[2]);
    trace(4,"pos=%.9f %.9f %.3f\n",pos[0]*R2D,pos[1]*R2D,pos[2]);
    for (i=0;i<n;i++) {
        trace(4,"sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
              obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
              azel[1+i*2]*R2D);
    }
    trace(4,"y=\n"); tracemat(4,y,nf*2,n,13,3);
    
    return 1;
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* DD (double-differenced) measurement error covariance ----------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;
    
    trace(3,"ddcov   : n=%d\n",n);
    
    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        
        // i==j时为Ri[k+i]+Rj[k+i]，i!=j时为Ri[k+i]
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,6);
}
/* baseline length constraint ------------------------------------------------*/
static int constbl(rtk_t *rtk, const double *x, const double *P, double *v,
                   double *H, double *Ri, double *Rj, int index)
{
    const double thres=0.1; /* threshold for nonliearity (v.2.3.0) */
    double xb[3],b[3],bb,var=0.0;
    int i;
     
    trace(3,"constbl : \n");
    
    /* no constraint */
    if (rtk->opt.baseline[0]<=0.0) return 0;
    
    /* time-adjusted baseline vector and length */
    for (i=0;i<3;i++) {
        xb[i]=rtk->rb[i];
        b[i]=x[i]-xb[i];
    }
    bb=norm(b,3);
    
    /* approximate variance of solution */
    if (P) {
        for (i=0;i<3;i++) var+=P[i+i*rtk->nx];
        var/=3.0;
    }
    /* check nonlinearity */
    if (var>SQR(thres*bb)) {
        trace(3,"constbl : equation nonlinear (bb=%.3f var=%.3f)\n",bb,var);
        return 0;
    }
    /* constraint to baseline length */
    v[index]=rtk->opt.baseline[0]-bb;
    if (H) {
        for (i=0;i<3;i++) H[i+index*rtk->nx]=b[i]/bb;
    }
    Ri[index]=0.0;
    Rj[index]=SQR(rtk->opt.baseline[1]);
    
    trace(4,"baseline len   v=%13.3f R=%8.6f %8.6f\n",v[index],Ri[index],Rj[index]);
    
    return 1;
}
/* precise tropspheric model -------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, int r,
                       const double *azel, const prcopt_t *opt, const double *x,
                       double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=IT(r,opt);
    
    // 调用 tropmapf() 计算干湿延迟投影系数
    /* wet mapping function */
    tropmapf(time,pos,azel,&m_w);   
    
    if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[i+1]+grad_e*x[i+2];
        dtdx[1]=grad_n*x[i];
        dtdx[2]=grad_e*x[i];
    }
    else dtdx[1]=dtdx[2]=0.0;
    dtdx[0]=m_w;
    return m_w*x[i];
}
/* test satellite system (m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN) ---------*/
static int test_sys(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
        case SYS_QZS: return m==4;
        case SYS_IRN: return m==5;
    }
    return 0;
}

// 计算观测方程系数矩阵 H、双差残差矩阵 v、双差观测值协方差矩阵 R
/* DD (double-differenced) phase/code residuals ------------------------------
 * args:rtk_t    *rtk      IO  rtk控制结构体
 *      nav_t    *nav      I   导航数据
 *      double   dt        I   接收机和基站的时间差
 *      double   *x        IO  状态变量
 *      double   *P        IO  状态变量的误差协方差阵
 *      int      sat       I   接收机和基站共同观测的卫星号列表
 *      double   *y        IO  相位/码残差
 *      double   *e        IO  观测矢量 (ecef)
 *      double   *azel     O   方位角和俯仰角 (rad)
 *      int      *iu       I   接收机和基站共同观测的卫星在接收机观测值中的index值列表
 *      int      *ir       I   接收机和基站共同观测的卫星在基站观测值中的index值列表
 *      int      ns        I   接收机和基站共同观测的卫星个数
 *      double   *v        O   实际观测量与预测观测量的残差
 *      double   *H        O   观测矩阵
 *      double   *R        O   测量误差的协方差
 *      int      *vflg     O   数据有效标志
 * return : (>0:ok,0:error)
 ------------------------------------------------------------------------------*/
static int ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                 const double *P, const int *sat, double *y, double *e,
                 double *azel, double *freq, const int *iu, const int *ir,
                 int ns, double *v, double *H, double *R, int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,freqi,freqj,*Hi=NULL;
    int i,j,k,m,f,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=NF(opt);
    
    trace(3,"ddres   : dt=%.1f nx=%d ns=%d\n",dt,rtk->nx,ns);
    
    bl=baseline(x,rtk->rb,dr);                  // 计算基线长度 bl，基线向量 dr
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);   // 基准站、流动站 XYZ-LLH
    
    Ri=mat(ns*nf*2+2,1); Rj=mat(ns*nf*2+2,1); im=mat(ns,1);
    tropu=mat(ns,1); tropr=mat(ns,1); dtdxu=mat(ns,3); dtdxr=mat(ns,3);
    
    // 所有伪距残差、载波相位残差重置为 0
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;  
    }
    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        // 若电离层模式 >=IONOOPT_EST，调用 ionmapf() 计算电离层延迟因子
        // 基准站、流动站，im[i]取基准站、流动站的平均
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(ionmapf(posu,azel+iu[i]*2)+ionmapf(posr,azel+ir[i]*2))/2.0;
        }
        // 若对流层模式 >=TROPOPT_EST，调用 prectrop() 计算对流层延迟因子
        // 流动站tropu[i]、基准站tropr[i]
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }

    // 遍历不同系统
    for (m=0;m<6;m++) /* m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN */
    
    // 遍历每一个频率，计算双差残差
    // 若为 PMODE_DGPS 伪距双差，需要限制遍历范围
    // 载波相位在 0-nf，伪距在 nf-2nf 为，因此伪距差分定位从 nf 开始
    for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) {
        // 寻找仰角最高的卫星作为参考卫星
        /* search reference satellite with highest elevation */
        for (i=-1,j=0;j<ns;j++) {
            sysi=rtk->ssat[sat[j]-1].sys;
            if (!test_sys(sysi,m)) continue;
            if (!validobs(iu[j],ir[j],f,nf,y)) continue;
            if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
        }
        if (i<0) continue;  // 选取失败则 continue
        
        // 遍历所有卫星，计算双差
        /* make DD (double difference) */
        for (j=0;j<ns;j++) {
            if (i==j) continue;
            sysi=rtk->ssat[sat[i]-1].sys;
            sysj=rtk->ssat[sat[j]-1].sys;
            freqi=freq[f%nf+iu[i]*nf];
            freqj=freq[f%nf+iu[j]*nf];
            if (!test_sys(sysj,m)) continue;
            if (!validobs(iu[j],ir[j],f,nf,y)) continue;
            
            // 用传入的非差相位/码残差 y 计算双差残差 v，并计算对应的 H
            if (H) {
                // 将 H+nv*rtk->nx 的地址赋给 H，此时更改 Hi 中的值即更改 H 矩阵中的值
                Hi=H+nv*rtk->nx;
                for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
            }

            //双差残差计算：【参考星(移动站) - 参考星(基准站)】- 【非参考星(移动站) - 非参考星(基准站)】
            /* DD residual */
            v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-
                  (y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);
            
            // 移动站位置偏导：移动站非参考星视线向量 - 移动站参考星视线向量
            /* partial derivatives by rover position */
            if (H) {
                for (k=0;k<3;k++) {
                    Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];   // 每一行前 3 列为观测向量差
                }
            }

            // 若要估计电离层参数，模式 IONOOPT_EST，用电离层延迟因子修正 v 和 H
            // 注意伪距和载波相位的电离层延迟大小相同，但是符号相反
            /* DD ionospheric delay term */
            if (opt->ionoopt==IONOOPT_EST) {
                didxi=(f<nf?-1.0:1.0)*im[i]*SQR(FREQ1/freqi);   // 载波相位为负
                didxj=(f<nf?-1.0:1.0)*im[j]*SQR(FREQ1/freqj);   // 伪距为正
                v[nv]-=didxi*x[II(sat[i],opt)]-didxj*x[II(sat[j],opt)];
                if (H) {
                    Hi[II(sat[i],opt)]= didxi;
                    Hi[II(sat[j],opt)]=-didxj;
                }
            }

            // 若要估计对流层参数，模式 TROPOPT_EST，用对流层延迟因子修正 v 和 H
            /* DD tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                    if (!H) continue;
                    Hi[IT(0,opt)+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
                    Hi[IT(1,opt)+k]=-(dtdxr[k+i*3]-dtdxr[k+j*3]);
                }
            }

            // 用相位偏移修正 v 和 H
            // 如果不是无电离层组合，从载波相位双差残差中扣除双差模糊度部分（即phase-bias），
            // 并对 H 矩阵中和模糊度相关的部分进行赋值
            /* DD phase-bias term */
            if (f<nf) {
                if (opt->ionoopt!=IONOOPT_IFLC) {
                    v[nv]-=CLIGHT/freqi*x[IB(sat[i],f,opt)]-
                           CLIGHT/freqj*x[IB(sat[j],f,opt)];
                    if (H) {
                        Hi[IB(sat[i],f,opt)]= CLIGHT/freqi;
                        Hi[IB(sat[j],f,opt)]=-CLIGHT/freqj;
                    }
                }
                else {
                    v[nv]-=x[IB(sat[i],f,opt)]-x[IB(sat[j],f,opt)];
                    if (H) {
                        Hi[IB(sat[i],f,opt)]= 1.0;
                        Hi[IB(sat[j],f,opt)]=-1.0;
                    }
                }
            }
            // 将双差载波，伪距残差，输出到相应的结构体中
            if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=v[nv];
            else      rtk->ssat[sat[j]-1].resp[f-nf]=v[nv];
            
            // 根据选项 maxinno 的值检测是否要排除此观测数据
            /* test innovation */
            if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                if (f<nf) {
                    rtk->ssat[sat[i]-1].rejc[f]++;
                    rtk->ssat[sat[j]-1].rejc[f]++;
                }
                errmsg(rtk,"outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                       sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                continue;
            }

            // 计算单差的测量误差协方差 Ri、i 为参考星，j 为非参考星
            /* SD (single-differenced) measurement error variances */
            Ri[nv]=varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt);
            Rj[nv]=varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt);
            
            // 设置数据有效标志，首先是载波，然后是伪距；根据 f<nf 判断
            /* set valid data flags */
            if (opt->mode>PMODE_DGPS) {
                if (f<nf) rtk->ssat[sat[i]-1].vsat[f]=rtk->ssat[sat[j]-1].vsat[f]=1;
            }
            else {
                rtk->ssat[sat[i]-1].vsat[f-nf]=rtk->ssat[sat[j]-1].vsat[f-nf]=1;
            }
            trace(4,"sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",sat[i],
                  sat[j],f<nf?"L":"P",f%nf+1,v[nv],Ri[nv],Rj[nv]);
            
            vflg[nv++]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
            nb[b]++;
        }
        b++;
    }
    /* end of system loop */
    
    // 如果是动基线的模式，增加基线长度约束
    /* baseline length constraint for moving baseline */
    if (opt->mode==PMODE_MOVEB&&constbl(rtk,x,P,v,H,Ri,Rj,nv)) {
        vflg[nv++]=3<<4;
        nb[b++]++;
    }
    if (H) {trace(5,"H=\n"); tracemat(5,H,rtk->nx,nv,7,4);}
    
    // 调用 ddcov()，计算载波相位/伪距双差量测噪声协方差阵 R
    /* DD measurement error covariance */
    ddcov(nb,b,Ri,Rj,nv,R);
    
    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr); free(dtdxu); free(dtdxr);
    
    return nv;
}
/* time-interpolation of residuals (for post-processing) ---------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2],freq[MAXOBS*NFREQ];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);
    
    trace(3,"intpres : n=%d tt=%.1f\n",n,tt);
    
    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;
    
    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);
    
    if (!zdres(1,obsb,nb,rs,dts,var,svh,nav,rtk->rb,opt,1,yb,e,azel,freq)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*nf*2,q=yb+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0; else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}
/* index for SD to DD transformation matrix D --------------------------------*/
static int ddidx(rtk_t *rtk, int *ix)
{
    int i,j,k,m,f,nb=0,na=rtk->na,nf=NF(&rtk->opt),nofix;
    
    trace(3,"ddidx   :\n");
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].fix[j]=0;
    }
    for (m=0;m<6;m++) { /* m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN */
        
        nofix=(m==1&&rtk->opt.glomodear==0)||(m==3&&rtk->opt.bdsmodear==0);
        
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {
            
            for (i=k;i<k+MAXSAT;i++) {
                if (rtk->x[i]==0.0||!test_sys(rtk->ssat[i-k].sys,m)||
                    !rtk->ssat[i-k].vsat[f]||!rtk->ssat[i-k].half[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    rtk->ssat[i-k].fix[f]=2; /* fix */
                    break;
                }
                else rtk->ssat[i-k].fix[f]=1;
            }
            for (j=k;j<k+MAXSAT;j++) {
                if (i==j||rtk->x[j]==0.0||!test_sys(rtk->ssat[j-k].sys,m)||
                    !rtk->ssat[j-k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[j-k].lock[f]>0&&!(rtk->ssat[j-k].slip[f]&2)&&
                    rtk->ssat[i-k].vsat[f]&&
                    rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    ix[nb*2  ]=i; /* state index of ref bias */
                    ix[nb*2+1]=j; /* state index of target bias */
                    nb++;
                    rtk->ssat[j-k].fix[f]=2; /* fix */
                }
                else rtk->ssat[j-k].fix[f]=1;
            }
        }
    }
    return nb;
}
/* restore SD (single-differenced) ambiguity ---------------------------------*/
static void restamb(rtk_t *rtk, const double *bias, int nb, double *xa)
{
    int i,n,m,f,index[MAXSAT],nv=0,nf=NF(&rtk->opt);
    
    trace(3,"restamb :\n");
    
    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];
    
    for (m=0;m<5;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            if (!test_sys(rtk->ssat[i].sys,m)||rtk->ssat[i].fix[f]!=2) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
        }
        if (n<2) continue;
        
        xa[index[0]]=rtk->x[index[0]];
        
        for (i=1;i<n;i++) {
            xa[index[i]]=xa[index[0]]-bias[nv++];
        }
    }
}

// 
/* hold integer ambiguity ----------------------------------------------------*/
static void holdamb(rtk_t *rtk, const double *xa)
{
    double *v,*H,*R;
    int i,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,nf=NF(&rtk->opt);
    
    trace(3,"holdamb :\n");
    
    v=mat(nb,1); H=zeros(nb,rtk->nx);
    
    // 循环遍历各个卫星，查找满足条件的卫星，并设置相应标志位 rtk->ssat[i].fix=3
    for (m=0;m<5;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            if (!test_sys(rtk->ssat[i].sys,m)||rtk->ssat[i].fix[f]!=2||
                rtk->ssat[i].azel[1]<rtk->opt.elmaskhold) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
            rtk->ssat[i].fix[f]=3; /* hold */
        }
        // 计算固定解双差和浮点解双差的差值，形成量测信息，并更新设计矩阵 H
        /* constraint to fixed ambiguity */
        for (i=1;i<n;i++) {
            v[nv]=(xa[index[0]]-xa[index[i]])-(rtk->x[index[0]]-rtk->x[index[i]]);
            
            H[index[0]+nv*rtk->nx]= 1.0;
            H[index[i]+nv*rtk->nx]=-1.0;
            nv++;
        }
    }
    // 若观测量数量有效，设置 R 阵，并调用 filter() 量测更新
    if (nv>0) {
        R=zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=VAR_HOLDAMB;
        
        /* update states with constraints */
        if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv))) {
            errmsg(rtk,"filter error (info=%d)\n",info);
        }
        free(R);
    }
    free(v); free(H);
}

// 根据 PRN 最小的卫星作为参考星，计算单差模糊度转双差模糊度的转换矩阵 D
// 输入双差模糊度及协方差矩阵，利用 lambda 算法固定模糊度。
// 根据矩阵反演公式，利用模糊度固定解调整其余状态参数浮点解，得到全部状态参数的固定解及协方差矩阵。
// 重塑单差模糊度，便于固定解有效性验证
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,nb,info,nx=rtk->nx,na=rtk->na;
    double *DP,*y,*b,*db,*Qb,*Qab,*QQ,s[2];
    int *ix;

    trace(3,"resamb_LAMBDA : nx=%d\n",nx);
    
    // 整周模糊度 ratio 赋初值0.0
    rtk->sol.ratio=0.0;
    
    // 检查配置中所设置的 ratio 值，如果该阈值小于1，则 return 0
    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }

    // 调用 ddidx() 函数，创建将卡尔曼状态量从单差转到双差的转换矩阵 D
    // 主要是将单差相位偏移状态量转换为双差相位偏移
    /* index of SD to DD transformation matrix D */
    ix=imat(nx,2);
    if ((nb=ddidx(rtk,ix))<=0) {
        errmsg(rtk,"no valid double-difference\n");
        free(ix);
        return 0;
    }

    // na 实际就是之前卡尔曼滤波中除了单差相位偏移之外的所有状态量个数（例如：位置+速度+加速度+电离层+对流层……）
    // nb 则是双差相位偏移的个数（即需要解算的整周模糊度个数）
    y=mat(nb,1); DP=mat(nb,nx-na); b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb);
    Qab=mat(na,nb); QQ=mat(na,nb);
    
    // 根据转移矩阵，求解双差整周模糊度以及协方差阵
    /* y=D*xc, Qb=D*Qc*D', Qab=Qac*D' */
    for (i=0;i<nb;i++) {
        y[i]=rtk->x[ix[i*2]]-rtk->x[ix[i*2+1]];
    }
    // 
    for (j=0;j<nx-na;j++) for (i=0;i<nb;i++) {
        DP[i+j*nb]=rtk->P[ix[i*2]+(na+j)*nx]-rtk->P[ix[i*2+1]+(na+j)*nx];
    }
    // 整周模糊度的协方差阵
    for (j=0;j<nb;j++) for (i=0;i<nb;i++) {
        Qb[i+j*nb]=DP[i+(ix[j*2]-na)*nb]-DP[i+(ix[j*2+1]-na)*nb];
    }
    // 整周和状态
    for (j=0;j<nb;j++) for (i=0;i<na;i++) {
        Qab[i+j*na]=rtk->P[i+ix[j*2]*nx]-rtk->P[i+ix[j*2+1]*nx];
    }

    // 调用 lambda() 函数计算双差整周模糊度最优解以及残差。估计结果存在 b、s 中
    /* LAMBDA/MLAMBDA ILS (integer least-square) estimation */
    if (!(info=lambda(nb,2,y,Qb,b,s))) {
        trace(4,"N(1)="); tracemat(4,b   ,1,nb,10,3);
        trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);
        
        // 计算 Ratio 值
        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

        // 如果 Ratio 值大于阈值，模糊度固定成功，求解固定解 rtk->xa 以及固定解的协方差 rtk->Pa
        /* validation by popular ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=opt->thresar[0]) {
            
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i];
                y[i]-=b[i];
            }
            if (!matinv(Qb,nb)) {
                matmul("NN",nb,1,nb, 1.0,Qb ,y,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db,1.0,rtk->xa);
                
                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa);
                
                trace(3,"resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
                
                // 调用 restamb() 重新存储单差模糊度
                /* restore SD ambiguity */
                restamb(rtk,bias,nb,xa);
            }
            else nb=0;
        }
        else { /* validation failed */
            errmsg(rtk,"ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb,s[1]/s[0],s[0],s[1]);
            nb=0;
        }
    }
    else {
        errmsg(rtk,"lambda error (info=%d)\n",info);
        nb=0;
    }
    free(ix);
    free(y); free(DP); free(b); free(db); free(Qb); free(Qab); free(QQ);
    
    return nb; /* number of ambiguities */
}
/* validation of solution ----------------------------------------------------*/
static int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                  int nv, double thres)
{
    double fact=thres*thres;
    int i,stat=1,sat1,sat2,type,freq;
    char *stype;
    
    trace(3,"valpos  : nv=%d thres=%.1f\n",nv,thres);
    
    /* post-fit residual test */
    for (i=0;i<nv;i++) {
        if (v[i]*v[i]<=fact*R[i+i*nv]) continue;
        sat1=(vflg[i]>>16)&0xFF;    // 参考卫星号
        sat2=(vflg[i]>> 8)&0xFF;    // 非参考卫星号
        type=(vflg[i]>> 4)&0xF;     // 种类：0载波 1伪距
        freq=vflg[i]&0xF;           // 第几个频率
        stype=type==0?"L":(type==1?"L":"C");
        errmsg(rtk,"large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
              sat1,sat2,stype,freq+1,v[i],SQRT(R[i+i*nv]));
    }
    return stat;
}

// 这个函数是 RTK 算法（实时动态相对定位）的核心函数
// 调用的 selsat、udstate、zdres、ddres、filter、resamb_LAMBDA 等函数都很复杂
/* relative positioning ------------------------------------------------------
 * args:rtk_t    *rtk      IO  rtk控制结构体
 *      obsd_t   *obs      I   观测数据
 *      int      nu        I   接收机观测数据的数量
 *      int      nr        I   基站观测数据的数量
 *      nav_t    *nav      I   导航数据
 * return : (1:ok,0:error)
 ----------------------------------------------------------------------------*/
static int relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                  const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*freq,*v,*H,*R,*xp,*Pp,*xa,*bias,dt;
    int i,j,f,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=opt->ionoopt==IONOOPT_IFLC?1:opt->nf;
    
    trace(3,"relpos  : nx=%d nu=%d nr=%d\n",rtk->nx,nu,nr);
    
    dt=timediff(time,obs[nu].time);     //计算流动站，参考站时间差
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); y=mat(nf*2,n); e=mat(3,n);
    azel=zeros(2,n); freq=zeros(nf,n);
    
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j]=0;
        for (j=1;j<NFREQ;j++) rtk->ssat[i].snr [j]=0;
    }

    // 根据卫星星历计算当前历元下各卫星的位置、速度、钟差
    /* satellite positions/clocks */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);
    
    // 计算基准站的各卫星观测值的非差残差（观测值-计算值）及卫星的高度角、方位角、卫星矢量等
    /* UD (undifferenced) residuals for base station */
    if (!zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,var+nu,svh+nu,nav,rtk->rb,opt,1,
               y+nu*nf*2,e+nu*3,azel+nu*2,freq+nu*nf)) {
        errmsg(rtk,"initial base station position error\n");
        
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        free(freq);
        return 0;
    }

    // 后处理中，需要时，调用 intpres() 进行插值
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt=intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }

    // 调用 selsat() 选择基准站和流动站之间的同步观测卫星，进行 RTK 算法时只需要基线间的同步观测卫星
    // 返回共同观测的卫星个数，输出卫星号列表 sat、在接收机观测值中的 index 值列表 iu 和在基站观测值中的 index 值列表 ir
    /* select common satellites between rover and base-station */
    if ((ns=selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        errmsg(rtk,"no common satellite\n");
        
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        free(freq);
        return 0;
    }

    // 调用 udstate 更新状态值 rtk->x 及其误差协方差 rtk->P
    /* temporal update of states */
    udstate(rtk,obs,sat,iu,ir,ns,nav);
    
    // 初始化变量内存以及赋初值
    trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,4);
    
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx); xa=mat(rtk->nx,1);
    matcpy(xp,rtk->x,rtk->nx,1);
    
    ny=ns*nf*2+2;
    v=mat(ny,1); H=zeros(rtk->nx,ny); R=mat(ny,ny); bias=mat(rtk->nx,1);
    
    // 设置迭代次数 niter，动基线加 2 次
    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);
    
    // 迭代量测更新 niter 次
    for (i=0;i<niter;i++) {
        // 计算流动站的各位卫星观测值非差残差（观测值-计算值）及卫星的高度角、方位角、卫星矢量等
        /* UD (undifferenced) residuals for rover */
        if (!zdres(0,obs,nu,rs,dts,var,svh,nav,xp,opt,0,y,e,azel,freq)) {
            errmsg(rtk,"rover initial position error\n");
            stat=SOLQ_NONE;
            break;
        }

        // 根据上述计算的基准站和流动站的各卫星观测值非差残差，计算双差残差矩阵 v
        // 并且在这个函数内，根据流动站非差卫星矢量计算观测方程的双差系数矩阵 H、
        // 根据非差观测值误差计算观测噪声的双差协方差矩阵 R
        /* DD (double-differenced) residuals and partial derivatives */
        if ((nv=ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,freq,iu,ir,ns,v,H,R,
                      vflg))<1) {
            errmsg(rtk,"no double-differenced residual\n");
            stat=SOLQ_NONE;
            break;
        }

        // 进行卡尔曼滤波的第二步——计算增益矩阵并状态更新，利用卡尔曼滤波计算 RTK 浮点解
        // 在进入 filter 函数之后需要对状态参数矩阵 x、状态噪声协方差矩阵 p
        // 观测方程系数阵H进行矩阵缩小简化，去除未用到的卫星
        // 这里 H 矩阵之所以使用转置的形式，是因为RTKLIB内存储矩阵的规则是按列存储的。
        /* Kalman filter measurement update */
        matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
        if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {
            errmsg(rtk,"filter error (info=%d)\n",info);
            stat=SOLQ_NONE;
            break;
        }
        trace(4,"x(%d)=",i+1); tracemat(4,xp,1,NR(opt),13,4);
    }

    // 量测更新完成，再次调用 zdres() 和 ddres() 计算双差相位/码残差，调用 valpos() 进行验证
    // 若通过则更新 rtk->x 以及 rtk->P，并更新模糊度控制结构体
    if (stat!=SOLQ_NONE&&zdres(0,obs,nu,rs,dts,var,svh,nav,xp,opt,0,y,e,azel,
                               freq)) {
        // 利用浮点结果计算双差残差和量测噪声
        /* post-fit residuals for float solution */
        nv=ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,freq,iu,ir,ns,v,NULL,R,vflg);
        
        // 计算流动站非差、双差残差，进行浮点解有效性验证
        /* validation of float solution */
        if (valpos(rtk,v,R,vflg,nv,4.0)) {

            /* update state and covariance matrix */    // 存储浮点结果
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);

            /* update ambiguity control struct */       // 存模糊度相关信息，统计有效卫星数
            rtk->sol.ns=0;
            for (i=0;i<ns;i++) for (f=0;f<nf;f++) {
                if (!rtk->ssat[sat[i]-1].vsat[f]) continue;
                rtk->ssat[sat[i]-1].lock[f]++;
                rtk->ssat[sat[i]-1].outc[f]=0;
                if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
            }

            /* lack of valid satellites */      // 检验卫星数是否足够
            if (rtk->sol.ns<4) stat=SOLQ_NONE;
        }
        else stat=SOLQ_NONE;
    }

    // 调用 resamb_LAMBDA() 利用 lambda 算法固定模糊度计算 RTK 固定解
    /* resolve integer ambiguity by LAMBDA */
    if (stat!=SOLQ_NONE&&resamb_LAMBDA(rtk,bias,xa)>1) {
        // 模糊度解算成功，根据固定结果计算残差和协方差，并进行校验
        if (zdres(0,obs,nu,rs,dts,var,svh,nav,xa,opt,0,y,e,azel,freq)) {
            
            /* post-fit reisiduals for fixed solution */
            nv=ddres(rtk,nav,dt,xa,NULL,sat,y,e,azel,freq,iu,ir,ns,v,NULL,R,
                     vflg);
            // 计算流动站非差、双差残差，进行固定解有效性验证
            /* validation of fixed solution */
            if (valpos(rtk,v,R,vflg,nv,4.0)) {
                // 固定解验证有效，若为 fix and hold 模式，需要存模糊度信息
                /* hold integer ambiguity */
                if (++rtk->nfix>=rtk->opt.minfix&&
                    rtk->opt.modear==ARMODE_FIXHOLD) {
                    holdamb(rtk,xa);
                }
                stat=SOLQ_FIX;
            }
        }
    }
    // 保存结果状态，位置、速度、方差
    /* save solution status */
    // 固定了存固定解，固定解在rtk->xa、rtk->Pa；否则存浮点解 rtk->x、rtk->P
    if (stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
        
        if (rtk->opt.dynamics) { /* velocity and covariance */
            for (i=3;i<6;i++) {
                rtk->sol.rr[i]=rtk->xa[i];
                rtk->sol.qv[i-3]=(float)rtk->Pa[i+i*rtk->na];
            }
            rtk->sol.qv[3]=(float)rtk->Pa[4+3*rtk->na];
            rtk->sol.qv[4]=(float)rtk->Pa[5+4*rtk->na];
            rtk->sol.qv[5]=(float)rtk->Pa[5+3*rtk->na];
        }
    }
    else {  // 浮点解数据在 rtk->x、rtk->P
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[1+2*rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
        
        if (rtk->opt.dynamics) { /* velocity and covariance */
            for (i=3;i<6;i++) {
                rtk->sol.rr[i]=rtk->x[i];
                rtk->sol.qv[i-3]=(float)rtk->P[i+i*rtk->nx];
            }
            rtk->sol.qv[3]=(float)rtk->P[4+3*rtk->nx];
            rtk->sol.qv[4]=(float)rtk->P[5+4*rtk->nx];
            rtk->sol.qv[5]=(float)rtk->P[5+3*rtk->nx];
        }
        rtk->nfix=0;
    }

    // 存当前历元载波信息，供下次使用
    for (i=0;i<n;i++) for (j=0;j<nf;j++) {
        if (obs[i].L[j]==0.0) continue;
        rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][j]=obs[i].time;    // 先前历元载波相位时间
        rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][j]=obs[i].L[j];    // 先前历元载波相位观测值
    }

    // 存 SNR 信噪比信息
    for (i=0;i<ns;i++) for (j=0;j<nf;j++) {
        
        /* output snr of rover receiver */
        rtk->ssat[sat[i]-1].snr[j]=obs[iu[i]].SNR[j];
    }

    // 存卫星的模糊度固定及周跳信息
    for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
        if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]++;
    }
    free(rs); free(dts); free(var); free(y); free(e); free(azel); free(freq);
    free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);
    
    if (stat!=SOLQ_NONE) rtk->sol.stat=stat;
    
    return stat!=SOLQ_NONE;
}

// 初始化 rtk_t 结构体
/* initialize RTK control ------------------------------------------------------
* initialize RTK control struct
* args   : rtk_t    *rtk    IO  TKk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    int i;
    
    trace(3,"rtkinit :\n");
    
    rtk->sol=sol0;
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=opt->mode<=PMODE_FIXED?NX(opt):pppnx(opt);
    rtk->na=opt->mode<=PMODE_FIXED?NR(opt):pppnx(opt);
    rtk->tt=0.0;
    rtk->x=zeros(rtk->nx,1);
    rtk->P=zeros(rtk->nx,rtk->nx);
    rtk->xa=zeros(rtk->na,1);
    rtk->Pa=zeros(rtk->na,rtk->na);
    rtk->nfix=rtk->neb=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
        rtk->ssat[i]=ssat0;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
    rtk->opt=*opt;
}

// 释放 rtk_t 结构体，释放给矩阵开辟的空间
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
{
    trace(3,"rtkfree :\n");
    
    rtk->nx=rtk->na=0;
    free(rtk->x ); rtk->x =NULL;
    free(rtk->P ); rtk->P =NULL;
    free(rtk->xa); rtk->xa=NULL;
    free(rtk->Pa); rtk->Pa=NULL;
}

// 这个函数的功能主要是实现实时动态定位（单历元定位），在这个函数内可选择不同的定位模式，包括 SPP，PPP，RTK 等  
/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by 
* precise positioning
* args   : rtk_t *rtk       IO  RTK control/result struct                           RTK 控制结构体
*            rtk->sol       IO  solution                                            结果结构体
*                .time      O   solution time                                       解算时间
*                .rr[]      IO  rover position/velocity                             流动站位置，速度，fixed 模式要输入，单点定位模式作为输出
*                               (I:fixed mode,O:single mode)                        
*                .dtr[0]    O   receiver clock bias (s)                             接收机 GPS 钟差
*                .dtr[1-5]  O   receiver GLO/GAL/BDS/IRN/QZS-GPS time offset (s)    接收机 GLO/GAL/BDS/IRN/QZS-GPS 钟差 (ISB)
*                .Qr[]      O   rover position covarinace                           流动站坐标协方差
*                .stat      O   solution status (SOLQ_???)                          解的状态SOLQ_???
*                .ns        O   number of valid satellites                          有效卫星数
*                .age       O   age of differential (s)                             差分龄期
*                .ratio     O   ratio factor for ambiguity validation               模糊度固定 Ratio 值
*            rtk->rb[]      IO  base station position/velocity                      基准站位置速度，相对定位一般做输入，动基线模式做输出
*                               (I:relative mode,O:moving-base mode)                
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellites in use                   有效卫星数
*            rtk->tt        O   time difference between current and previous (s)    相邻历元间隔
*            rtk->x[]       IO  float states pre-filter and post-filter             浮点解
*            rtk->P[]       IO  float covariance pre-filter and post-filter         浮点解协方差
*            rtk->xa[]      O   fixed states after AR                               固定解 
*            rtk->Pa[]      O   fixed covariance after AR                           固定解协方差
*            rtk->ssat[s]   IO  satellite {s+1} status                              卫星结果状态控制结构体数组
*                .sys       O   system (SYS_???)                                    卫星导航系统 
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)            方位角 
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)            高度角
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)            
*                .resp [f]  O   freq(f+1) pseudorange residual (m)                  伪距残差 
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)                载波相位残差
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)            
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) cycle slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
extern int rtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    sol_t solb={{0}};
    gtime_t time;
    int i,nu,nr;
    char msg[128]="";
    
    trace(3,"rtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
    trace(4,"obs=\n"); traceobs(4,obs,n);
    
    // 设置 rtk 内基准站坐标，基准站坐标在 execses() 函数内已经计算了，速度设为 0.0
    // 这里将配置结构体 opt 内基准站的坐标赋值给解算结构体 rtk 内基准站的坐标
    /* set base staion position */
    if (opt->refpos<=POSOPT_RINEX&&opt->mode!=PMODE_SINGLE&&
        opt->mode!=PMODE_MOVEB) {
        for (i=0;i<6;i++) rtk->rb[i]=i<3?opt->rb[i]:0.0;
    }

    // 统计基准站、流动站观测值个数 nu、nr，可用于后面判断是否满足差分条件
    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;
    
    time=rtk->sol.time; /* previous epoch */
    
    // 利用观测值及星历计算流动站的 SPP 定位结果，作为 kalman 滤波的近似坐标。需要注意，
    // 如果由于流动站 SPP 定位结果坐标误差过大等原因导致的 SPP 无解，则不进行 rtk 运算，当前历元无解
    /* rover position by single point positioning */
    if (!pntpos(obs,nu,nav,&rtk->opt,&rtk->sol,NULL,rtk->ssat,msg)) {
        errmsg(rtk,"point pos error (%s)\n",msg);
        
        if (!rtk->opt.dynamics) {
            outsolstat(rtk);
            return 0;
        }
    }

    // 计算当前历元和上一历元时间差 rtk->tt，其中 rtk->sol.time 是当前历元时间，time 是上一历元时间
    if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);
    
    // 如果是单点定位模式直接输出刚刚 SPP 算的坐标
    /* single point positioning */
    if (opt->mode==PMODE_SINGLE) {
        outsolstat(rtk);
        return 1;
    }

    // 如果不是单点模式，抑制单点解的输出
    /* suppress output of single solution */
    if (!opt->outsingle) {
        rtk->sol.stat=SOLQ_NONE;
    }

    // 精密单点定位模式调用 pppos() 进行解算
    /* precise point positioning */
    if (opt->mode>=PMODE_PPP_KINEMA) {
        pppos(rtk,obs,nu,nav);
        outsolstat(rtk);
        return 1;
    }

    // 程序能运行到这，一定是相对定位模式了
    // 检查该历元流动站观测时间和基准站观测时间是否对应，若无基准站观测数据直接 return
    /* check number of data of base station and age of differential */
    if (nr==0) {
        errmsg(rtk,"no base station observation data for rtk\n");
        outsolstat(rtk);
        return 1;
    }

    //  动基线的基站坐标需要随时间同步变化，所以需要计算出变化速率,
    //  解释了前面除了单点定位，动基线也不参与基站解算，动基线在这里单独解算
    if (opt->mode==PMODE_MOVEB) { /*  moving baseline */
        
        // 调用 pntpos() 计算基准站位置
        /* estimate position/velocity of base station */
        if (!pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
            errmsg(rtk,"base station position error (%s)\n",msg);
            return 0;
        }

        // 计算差分龄期 rtk->sol.age
        rtk->sol.age=(float)timediff(rtk->sol.time,solb.time);
        
        if (fabs(rtk->sol.age)>TTOL_MOVEB) {
            errmsg(rtk,"time sync error for moving-base (age=%.1f)\n",rtk->sol.age);
            return 0;
        }
        for (i=0;i<6;i++) rtk->rb[i]=solb.rr[i];
        
        // 把 solb.rr 赋值给 rtk->rb
        /* time-synchronized position of base station */
        for (i=0;i<3;i++) rtk->rb[i]+=rtk->rb[i+3]*rtk->sol.age;
    }
    else {
        rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);
        
        if (fabs(rtk->sol.age)>opt->maxtdiff) {
            errmsg(rtk,"age of differential error (age=%.1f)\n",rtk->sol.age);
            outsolstat(rtk);
            return 1;
        }
    }

    // 上面的步骤只算了相对定位的差分龄期和动基线坐标，这里调用 relpos() 进行相位定位，并输出最终结果
    /* relative potitioning */
    relpos(rtk,obs,nu,nr,nav);
    outsolstat(rtk);
    
    return 1;
}
