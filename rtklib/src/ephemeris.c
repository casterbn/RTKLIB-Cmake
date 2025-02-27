/*------------------------------------------------------------------------------
* ephemeris.c : satellite ephemeris and clock functions
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] IS-GPS-200K, Navstar GPS Space Segment/Navigation User Interfaces,
*         May 6, 2019
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Version 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performance standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1.3, December, 2016
*     [8] Quasi-Zenith Satellite System Interface Specification Satellite
*         Positioning, Navigation and Timing Service (IS-QZSS-PNT-003), Cabinet
*         Office, November 5, 2018
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 3.0), China Satellite
*         Navigation office, February, 2019
*     [10] RTCM Standard 10403.3, Differential GNSS (Global Navigation
*         Satellite Systems) Services - version 3, October 7, 2016
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.1  moved from rtkcmn.c
*                           added api:
*                               eph2clk(),geph2clk(),seph2clk(),satantoff()
*                               satposs()
*                           changed api:
*                               eph2pos(),geph2pos(),satpos()
*                           deleted api:
*                               satposv(),satposiode()
*           2010/08/26 1.2  add ephemeris option EPHOPT_LEX
*           2010/09/09 1.3  fix problem when precise clock outage
*           2011/01/12 1.4  add api alm2pos()
*                           change api satpos(),satposs()
*                           enable valid unhealthy satellites and output status
*                           fix bug on exception by glonass ephem computation
*           2013/01/10 1.5  support beidou (compass)
*                           use newton's method to solve kepler eq.
*                           update ssr correction algorithm
*           2013/03/20 1.6  fix problem on ssr clock relativitic correction
*           2013/09/01 1.7  support negative pseudorange
*                           fix bug on variance in case of ura ssr = 63
*           2013/11/11 1.8  change constant MAXAGESSR 70.0 -> 90.0
*           2014/10/24 1.9  fix bug on return of var_uraeph() if ura<0||15<ura
*           2014/12/07 1.10 modify MAXDTOE for qzss,gal and bds
*                           test max number of iteration for Kepler
*           2015/08/26 1.11 update RTOL_ELPLER 1E-14 -> 1E-13
*                           set MAX_ITER_KEPLER for alm2pos()
*           2017/04/11 1.12 fix bug on max number of obs data in satposs()
*           2018/10/10 1.13 update reference [7]
*                           support ura value in var_uraeph() for galileo
*                           test eph->flag to recognize beidou geo
*                           add api satseleph() for ephemeris selection
*           2020/11/30 1.14 update references [1],[2],[8],[9] and [10]
*                           add API getseleph()
*                           rename API satseleph() as setseleph()
*                           support NavIC/IRNSS by API satpos() and satposs()
*                           support BDS C59-63 as GEO satellites in eph2pos()
*                           default selection of I/NAV for Galileo ephemeris
*                           no support EPHOPT_LEX by API satpos() and satposs()
*                           unselect Galileo ephemeris with AOD<=0 in seleph()
*                           fix bug on clock iteration in eph2clk(), geph2clk()
*                           fix bug on clock reference time in satpos_ssr()
*                           fix bug on wrong value with ura=15 in var_ura()
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

// 用 URA 用户测距精度标定方差
/* variance by ura ephemeris -------------------------------------------------*/
static double var_uraeph(int sys, int ura)
{
    const double ura_value[]={   
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||14<ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
static double var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}

// 用 alm_t 中存的历书求卫星位置、钟差
/* almanac to satellite position and clock bias --------------------------------
* compute satellite position and clock bias with almanac (gps, galileo, qzss)
* args   : gtime_t time     I   信号发射时间 (gpst)
*          alm_t *alm       I   历书参数结构体
*          double *rs       O   卫星位置 (ecef) {x,y,z} (m)
*          double *dts      O   卫星钟差 (s)
* return : none
* notes  : see ref [1],[7],[8]
*-----------------------------------------------------------------------------*/
extern void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,x,y,sinO,cosO,cosi,mu;
    int n;
    
    trace(4,"alm2pos : time=%s sat=%2d\n",time_str(time,3),alm->sat);
    
    tk=timediff(time,alm->toa);
    
    if (alm->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=0.0;
        return;
    }
    mu=satsys(alm->sat,NULL)==SYS_GAL?MU_GAL:MU_GPS;
    
    M=alm->M0+sqrt(mu/(alm->A*alm->A*alm->A))*tk;
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-alm->e*sin(E)-M)/(1.0-alm->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        trace(2,"alm2pos: kepler iteration overflow sat=%2d\n",alm->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);
    u=atan2(sqrt(1.0-alm->e*alm->e)*sinE,cosE-alm->e)+alm->omg;
    r=alm->A*(1.0-alm->e*cosE);
    i=alm->i0;
    O=alm->OMG0+(alm->OMGd-OMGE)*tk-OMGE*alm->toas;
    x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);
    rs[0]=x*cosO-y*cosi*sinO;
    rs[1]=x*sinO+y*cosi*cosO;
    rs[2]=y*sin(i);
    *dts=alm->f0+alm->f1*tk;
}

// 根据信号发射时间和广播星历，计算卫星钟差，不考虑相对论效应和 TGD
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
    double t,ts;
    int i;
    
    trace(4,"eph2clk : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
    // 计算与星历参考时间的偏差 dt = t-toc
    t=ts=timediff(time,eph->toc);
    
    // 利用二项式校正计算出卫星钟差，从 dt 中减去这部分，然后再进行一次上述操作，得到最终的 dt
    for (i=0;i<2;i++) {
        t=ts-(eph->f0+eph->f1*t+eph->f2*t*t);   // (E.4.16)
    }

    // 使用二项式校正得到最终的钟差
    return eph->f0+eph->f1*t+eph->f2*t*t;
}

// 根据广播星历计算出算信号发射时刻卫星的位置和钟差
/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
                    double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    
    trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
    // 通过卫星轨道半长轴 A 判断星历是否有效，无效则返回
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    // 计算规化时间 tk (E.4.2)
    tk=timediff(time,eph->toe);
    
    // 根据不同卫星系统设置相应的地球引力常数 mu 和 地球自转角速度 omge
    switch ((sys=satsys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }

    // 计算平近点角 M (E.4.3)
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    
    // 用牛顿迭代法来计算偏近点角 E。参考 RTKLIB manual P145 (E.4.19) (E.4.4)
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        trace(2,"eph2pos: kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);
    
    trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);
    
    // 计算摄动改正后的升交点角距 u 卫星矢径长度 r 轨道倾角 i
    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;     // (E.4.5) (E.4.6) (E.4.10)
    r=eph->A*(1.0-eph->e*cosE);             // (E.4.11)
    i=eph->i0+eph->idot*tk;                 // (E.4.12)
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);     
    u+=eph->cus*sin2u+eph->cuc*cos2u;       // (E.4.7)
    r+=eph->crs*sin2u+eph->crc*cos2u;       // (E.4.8)
    i+=eph->cis*sin2u+eph->cic*cos2u;       // (E.4.9)
    x=r*cos(u); y=r*sin(u); cosi=cos(i);
    
    // 北斗的 MEO、IGSO 卫星计算方法与 GPS, Galileo and QZSS 相同，只是一些参数不同
    // GEO 卫星的 O 和最后位置的计算稍有不同 
    /* beidou geo satellite */
    if (sys==SYS_CMP&&(prn<=5||prn>=59)) { /* ref [9] table 4-1 */
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;    // (E.4.29)
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5; // ECEF位置(E.4.30)
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes; // 计算升交点赤经O (E.4.13)
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;       // 计算卫星ECEF位置存入 rs 中 (E.4.14)
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);             // (E.4.15)
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;  // 利用三个二项式模型系数 af0、af1、af2 计算卫星钟差
    
    /* relativity correction */         // 相对论效应改正卫星钟差
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
    
    /* position and clock error variance */
    *var=var_uraeph(sys,eph->sva);      //用 URA 值来标定方差
}

// GLONASS 轨道微分方程，xdot[3]:dvx/dt，xdot[4]:dvy/dt
/* glonass orbit differential equations --------------------------------------
 * args:x[6]    I   上一次计算时卫星位置速度
 *      xdot[6] O   差分方程
 *      acc[3]  I   加速度
 ----------------------------------------------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,
        r2=dot(x,x,3),          // r平方
        r3=r2*sqrt(r2),         // r三次方
        omg2=SQR(OMGE_GLO);     // omg平方
    
    if (r2<=0.0) {      // 计算出错
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];           // (E.4.22)
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];     // (E.4.23)
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];     // (E.4.24)
    xdot[5]=(c-2.0*a)*x[2]+acc[2];                      // (E.4.25)
}

// 数值积分计算 GLONASS 位置、速度
/* glonass position and velocity by numerical integration --------------------*/
static void glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;
    
    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
    double t,ts;
    int i;
    
    trace(4,"geph2clk: time=%s sat=%2d\n",time_str(time,3),geph->sat);
    
    t=ts=timediff(time,geph->toe);
    
    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);    // (E.4.26)
    }
    return -geph->taun+geph->gamn*t;
}

/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var)
{
    double t,tt,x[6];
    int i;
    
    trace(4,"geph2pos: time=%s sat=%2d\n",time_str(time,3),geph->sat);
    
    t=timediff(time,geph->toe);
    
    *dts=-geph->taun+geph->gamn*t;      // 计算钟差dts(E.4.26)
    
    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }

    // 步长 TSTEP：60s
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];
    
    // GLONASS 卫星的方差直接定为 5*5
    *var=SQR(ERREPH_GLO);
}

/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern double seph2clk(gtime_t time, const seph_t *seph)
{
    double t;
    int i;
    
    trace(4,"seph2clk: time=%s sat=%2d\n",time_str(time,3),seph->sat);
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<2;i++) {             // 由参数 af0/af1，钟差校正
        t-=seph->af0+seph->af1*t;   // (E.4.33)
    }
    return seph->af0+seph->af1*t;
}
/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args   : gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var)
{
    double t;
    int i;
    
    trace(4,"seph2pos: time=%s sat=%2d\n",time_str(time,3),seph->sat);
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<3;i++) {             // 由位置，速度，加速度计算位置 (E.4.32)
        rs[i]=seph->pos[i]+seph->vel[i]*t+seph->acc[i]*t*t/2.0;
    }
    *dts=seph->af0+seph->af1*t;     // 钟差 (E.4.33)
    
    *var=var_uraeph(SYS_SBS,seph->sva);
}

// 传入 sattle number，时间或 IDOE，如果传入 IDOE>=0，按 IDOE 找星历数据，否则取最接近时间的星历
/* select ephememeris --------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax,tmin;
    int i,j=-1,sys,sel;
    
    trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    // 根据传入的 sattle number，调用 satsys() 判断卫星系统，赋值 tmax，tmin，sel
    sys=satsys(sat,NULL);
    switch (sys) {
        case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=eph_sel[0]; break;
        case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=eph_sel[2]; break;
        case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=eph_sel[3]; break;
        case SYS_CMP: tmax=MAXDTOE_CMP+1.0; sel=eph_sel[4]; break;
        case SYS_IRN: tmax=MAXDTOE_IRN+1.0; sel=eph_sel[5]; break;
        default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;
    
    // 遍历 nav->eph[]
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;                 // eph[i] 不是需要的卫星，进入下一次循环
        if (iode>=0&&nav->eph[i].iode!=iode) continue;      // 如果传入了 idoe 时间不符合，也进行下一次循环
        
        // 若是伽利略卫星还要判断是 I/NAV、F/NAV
        if (sys==SYS_GAL) {
            sel=getseleph(SYS_GAL);
            if (sel==0&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==1&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
            if (timediff(nav->eph[i].toe,time)>=0.0) continue; /* AOD<=0 */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;    // 时间差过大，也进行下一次循化
        if (iode>=0) return nav->eph+i;                                 // 传入IDOE>=0，直接返回符合条件的星历
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */           // 存下时间差最小的星历下标和时间差
    }
    if (iode>=0||j<0) {
        trace(3,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",
              time_str(time,0),sat,iode);
        return NULL;
    }
    return nav->eph+j;
}
/* select glonass ephememeris ------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;
    
    trace(4,"selgeph : time=%s sat=%2d iode=%2d\n",time_str(time,3),sat,iode);
    
    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;                            // 卫星不同，进行下一次循环
        if (iode>=0&&nav->geph[i].iode!=iode) continue;                 // 传入 IDOE>0，IDOE 不同，进行下一次循环
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;   
        if (iode>=0) return nav->geph+i;                                // 传入IDOE>=0，直接返回符合条件的星历
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */           // 存下时间差最小的星历下标和时间差
    }
    if (iode>=0||j<0) {
        trace(3,"no glonass ephemeris  : %s sat=%2d iode=%2d\n",time_str(time,0),
              sat,iode);
        return NULL;
    }
    return nav->geph+j;
}
/* select sbas ephememeris ---------------------------------------------------*/
static seph_t *selseph(gtime_t time, int sat, const nav_t *nav)
{
    double t,tmax=MAXDTOE_SBS,tmin=tmax+1.0;
    int i,j=-1;
    
    trace(4,"selseph : time=%s sat=%2d\n",time_str(time,3),sat);
    
    for (i=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=sat) continue;
        if ((t=fabs(timediff(nav->seph[i].t0,time)))>tmax) continue;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (j<0) {
        trace(3,"no sbas ephemeris     : %s sat=%2d\n",time_str(time,0),sat);
        return NULL;
    }
    return nav->seph+j;
}

// 通过广播星历来确定卫星钟差
/* satellite clock with broadcast ephemeris ----------------------------------
 * arge:gtime_t  time      I   信号发射时刻
 *      gtime_t  teph      I   用于选择星历的时刻 (gpst)
 *      int      sat       I   卫星号 (1-MAXSAT)
 *      nav_t    *nav      I   导航数据
 *      double   *dts      O   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 * return:(1:ok,0:error)
----------------------------------------------------------------------------*/
static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    int sys;
    
    trace(4,"ephclk  : time=%s sat=%2d\n",time_str(time,3),sat);

    //调用 satsys() 函数，根据卫星编号确定该卫星所属的导航系统和该卫星在该系统中的 PRN编号
    sys=satsys(sat,NULL);
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        // 调用 seleph() 函数来选择最接近 teph 的那个星
        if (!(eph=seleph(teph,sat,-1,nav))) return 0;
        // 调用 eph2clk() 函数，通过广播星历和信号发射时间计算出卫星钟差
        *dts=eph2clk(time,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,geph);
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        *dts=seph2clk(time,seph);
    }
    else return 0;
    
    return 1;
}

// 根据卫星系统，先调用对应的星历选择函数，seleph()、selgeph()、seph2pos()；
// 再调用对应的解算函数，eph2pos()、geph2pos()、seph2pos()计算位置、钟差
// 增加一个极短的时间tt，再调用对应的解算函数计算位置、钟差两次的位置、钟差相减再除以tt，得速度、钟漂
/* satellite position and clock by broadcast ephemeris -----------------------
 * arge:gtime_t  time      I   信号发射时刻
 *      gtime_t  teph      I   用于选择星历的时刻 (gpst)
 *      int      sat       I   卫星号 (1-MAXSAT)
 *      nav_t    *nav      I   导航数据
 *      int      iode      I   星历数据期号
 *      double   *rs       O   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
 *      double   *dts      O   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
 *      double   *var      O   卫星位置和钟差的协方差 (m^2)
 *      int      *svh      O   卫星健康标志 (-1:correction not available)
 * return:(1:ok,0:error)
----------------------------------------------------------------------------*/
static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  int iode, double *rs, double *dts, double *var, int *svh)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;
    
    trace(4,"ephpos  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    // 调用 satsys() 函数，确定该卫星所属的导航系统
    sys=satsys(sat,NULL);
    
    *svh=-1;
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        // 调用 seleph() 函数来选择广播星历
        if (!(eph=seleph(teph,sat,iode,nav))) return 0;
        // 根据选中的广播星历，调用 eph2pos() 函数来计算信号发射时刻卫星的位置、钟差和相应结果的误差
        eph2pos(time,eph,rs,dts,var);
        time=timeadd(time,tt);
        eph2pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        seph2pos(time,seph,rs,dts,var);
        time=timeadd(time,tt);
        seph2pos(time,seph,rst,dtst,var);
        *svh=seph->svh;
    }
    else return 0;
    
    // 在信号发射时刻的基础上给定一个微小的时间间隔，再次计算新时刻的 P、V、C
    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;    // 卫星速度 rs[i+3]
    dts[1]=(dtst[0]-dts[0])/tt;                     // 钟漂 dts[1]
    
    return 1;
}
/* satellite position and clock with sbas correction -------------------------*/
static int satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                        double *rs, double *dts, double *var, int *svh)
{
    const sbssatp_t *sbs;
    int i;
    
    trace(4,"satpos_sbas: time=%s sat=%2d\n",time_str(time,3),sat);
    
    // 寻找要计算卫星对应的SBAS信息
    /* search sbas satellite correciton */      
    for (i=0;i<nav->sbssat.nsat;i++) {
        sbs=nav->sbssat.sat+i;
        if (sbs->sat==sat) break;
    }

    // 找不到，则 svh 写 -1，只调用 ephpos() 用广播星历计算卫星位置
    if (i>=nav->sbssat.nsat) {                  
        trace(2,"no sbas correction for orbit: %s sat=%2d\n",time_str(time,0),sat);
        ephpos(time,teph,sat,nav,-1,rs,dts,var,svh);
        *svh=-1;
        return 0;
    }

    // 用广播星历计算卫星位置
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time,teph,sat,nav,sbs->lcorr.iode,rs,dts,var,svh)) return 0;
    
    /* sbas satellite correction (long term and fast) */
    if (sbssatcorr(time,sat,nav,rs,dts,var)) return 1;
    *svh=-1;
    return 0;
}
/* satellite position and clock with ssr correction --------------------------*/
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                      int opt, double *rs, double *dts, double *var, int *svh)
{
    const ssr_t *ssr;
    eph_t *eph;
    double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk;
    int i,sys;
    
    trace(4,"satpos_ssr: time=%s sat=%2d\n",time_str(time,3),sat);
    
    ssr=nav->ssr+sat-1;
    
    if (!ssr->t0[0].time) {
        trace(2,"no ssr orbit correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (!ssr->t0[1].time) {
        trace(2,"no ssr clock correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* inconsistency between orbit and clock correction */
    if (ssr->iod[0]!=ssr->iod[1]) {
        trace(2,"inconsist ssr correction: %s sat=%2d iod=%d %d\n",
              time_str(time,0),sat,ssr->iod[0],ssr->iod[1]);
        *svh=-1;
        return 0;
    }
    t1=timediff(time,ssr->t0[0]);
    t2=timediff(time,ssr->t0[1]);
    t3=timediff(time,ssr->t0[2]);
    
    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1)>MAXAGESSR||fabs(t2)>MAXAGESSR) {
        trace(2,"age of ssr error: %s sat=%2d t=%.0f %.0f\n",time_str(time,0),
              sat,t1,t2);
        *svh=-1;
        return 0;
    }
    if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
    if (ssr->udi[1]>=1.0) t2-=ssr->udi[1]/2.0;
    
    for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
    dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;
    
    /* ssr highrate clock correction (ref [4]) */
    if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) {
        dclk+=ssr->hrclk;
    }
    if (norm(deph,3)>MAXECORSSR||fabs(dclk)>MAXCCORSSR) {
        trace(3,"invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
              time_str(time,0),norm(deph,3),dclk);
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time,teph,sat,nav,ssr->iode,rs,dts,var,svh)) return 0;
    
    /* satellite clock for gps, galileo and qzss */
    sys=satsys(sat,NULL);
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=seleph(teph,sat,ssr->iode,nav))) return 0;
        
        /* satellite clock by clock parameters */
        tk=timediff(time,eph->toc);
        dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
        dts[1]=eph->f1+2.0*eph->f2*tk;
        
        /* relativity correction */
        dts[0]-=2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
    }
    /* radial-along-cross directions in ecef */
    if (!normv3(rs+3,ea)) return 0;
    cross3(rs,rs+3,rc);
    if (!normv3(rc,ec)) {
        *svh=-1;
        return 0;
    }
    cross3(ea,ec,er);
    
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rs,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i]+=-(er[i]*deph[0]+ea[i]*deph[1]+ec[i]*deph[2])+dant[i];
    }
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    dts[0]+=dclk/CLIGHT;
    
    /* variance by ssr ura */
    *var=var_urassr(ssr->ura);
    
    trace(5,"satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
          time_str(time,2),sat,deph[0],deph[1],deph[2],er[0],er[1],er[2],dclk,*var);
    
    return 1;
}

// 计算信号发射时刻卫星的 P(ecef,m)、V(ecef,m/s)、C((s|s/s))
/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock                                
* args   : gtime_t time     I   time (gpst)                                     卫星信号发射时间
*          gtime_t teph     I   time to select ephemeris (gpst)                 用于选择星历的时刻
*          int    sat       I   satellite number                                卫星号
*          nav_t  *nav      I   navigation data                                 NAV 导航数据
*          int    ephopt    I   ephemeris option (EPHOPT_???)                   星历选项 (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)                卫星位置和速度，长度为 6*n
*                               {x,y,z,vx,vy,vz} (m|m/s)                        
*          double *dts      O   sat clock {bias,drift} (s|s/s)                  卫星钟差，长度为2*n， {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)     卫星位置和钟差的协方差 (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)   卫星健康标志 (-1:correction not available)
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                  const nav_t *nav, double *rs, double *dts, double *var,
                  int *svh)
{
    trace(4,"satpos  : time=%s sat=%2d ephopt=%d\n",time_str(time,3),sat,ephopt);
    
    *svh=0;
    
    // 根据星历选项调用对应的解算函数
    switch (ephopt) {
        case EPHOPT_BRDC  : return ephpos     (time,teph,sat,nav,-1,rs,dts,var,svh);    // 广播星历
        case EPHOPT_SBAS  : return satpos_sbas(time,teph,sat,nav,   rs,dts,var,svh);    // SBAS
        case EPHOPT_SSRAPC: return satpos_ssr (time,teph,sat,nav, 0,rs,dts,var,svh);    // 参考天线相位中心
        case EPHOPT_SSRCOM: return satpos_ssr (time,teph,sat,nav, 1,rs,dts,var,svh);    // 参考质心，还需要天线相位中心改正
        case EPHOPT_PREC  :
            if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;             // 精密星历
    }
    *svh=-1;
    return 0;
}

// 按照所观测到的卫星顺序计算出每颗卫星的位置、速度、{钟差、频漂}
/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : gtime_t teph     I   time to select ephemeris (gpst)                 (gpst) 用于选择星历的时刻 (gpst)
*          obsd_t *obs      I   observation data                                OBS 观测数据
*          int    n         I   number of observation data                      OBS 数
*          nav_t  *nav      I   navigation data                                 NAV 导航电文
*          int    ephopt    I   ephemeris option (EPHOPT_???)                   星历选项 (EPHOPT_???)
*          double *rs       O   satellite positions and velocities (ecef)       卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
*          double *dts      O   satellite clocks                                卫星钟差，长度为2*n， {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variances (m^2)    卫星位置和钟差的协方差 (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)   卫星健康标志 (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
                    int ephopt, double *rs, double *dts, double *var, int *svh)
{
    gtime_t time[2*MAXOBS]={{0}};
    double dt,pr;
    int i,j;
    
    trace(3,"satposs : teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,ephopt);
    
    // 遍历每一条观测数据
    for (i=0;i<n&&i<2*MAXOBS;i++) {

        // 首先初始化，将对当前观测数据的 rs、dts、var和svh数组的元素置 0
        for (j=0;j<6;j++) rs [j+i*6]=0.0;
        for (j=0;j<2;j++) dts[j+i*2]=0.0;
        var[i]=0.0; svh[i]=0;
        
        // 通过判断某一频率下信号的伪距是否为 0，来得到此时所用的频率个数
        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) {
            trace(3,"no pseudorange %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
            continue;
        }

        // 用数据接收时刻减去伪距信号传播时间，得到卫星信号的发射时刻
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);
        
        // 调用 ephclk() 函数，由广播星历计算出当前观测卫星与 GPS 时间的钟差 dt
        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)) {
            trace(3,"no broadcast clock %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }
        // 信号发射时刻减去钟差 dt，得到 GPS 时间下的卫星信号发射时刻
        time[i]=timeadd(time[i],-dt);
        
        // 调用 satpos() 函数，计算信号发射时刻卫星的位置(ecef,m)、速度(ecef,m/s)、钟差((s|s/s))
        // 注意，这里计算出的钟差是考虑了相对论效应的了，只是还没有考虑 TGD
        /* satellite position and clock at transmission time */
        if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*6,dts+i*2,var+i,
                    svh+i)) {
            trace(3,"no ephemeris %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }

        // 如果没有精密钟差，则用广播星历的钟差替代，钟漂设为 0
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*2]==0.0) {
            if (!ephclk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
            dts[1+i*2]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
              time_str(time[i],6),obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],
              dts[i*2]*1E9,var[i],svh[i]);
    }
}
/* set selected satellite ephemeris --------------------------------------------
* Set selected satellite ephemeris for multiple ones like LNAV - CNAV, I/NAV -
* F/NAV. Call it before calling satpos(),satposs() to use unselected one.
* args   : int    sys       I   satellite system (SYS_???)
*          int    sel       I   selection of ephemeris
*                                 GPS,QZS : 0:LNAV ,1:CNAV  (default: LNAV)
*                                 GAL     : 0:I/NAV,1:F/NAV (default: I/NAV)
*                                 others  : undefined
* return : none
* notes  : default ephemeris selection for galileo is any.
*-----------------------------------------------------------------------------*/
extern void setseleph(int sys, int sel)
{
    switch (sys) {
        case SYS_GPS: eph_sel[0]=sel; break;
        case SYS_GLO: eph_sel[1]=sel; break;
        case SYS_GAL: eph_sel[2]=sel; break;
        case SYS_QZS: eph_sel[3]=sel; break;
        case SYS_CMP: eph_sel[4]=sel; break;
        case SYS_IRN: eph_sel[5]=sel; break;
        case SYS_SBS: eph_sel[6]=sel; break;
    }
}
/* get selected satellite ephemeris -------------------------------------------
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
extern int getseleph(int sys)
{
    switch (sys) {
        case SYS_GPS: return eph_sel[0];
        case SYS_GLO: return eph_sel[1];
        case SYS_GAL: return eph_sel[2];
        case SYS_QZS: return eph_sel[3];
        case SYS_CMP: return eph_sel[4];
        case SYS_IRN: return eph_sel[5];
        case SYS_SBS: return eph_sel[6];
    }
    return 0;
}
