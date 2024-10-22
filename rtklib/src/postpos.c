/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioing api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*           2015/03/23  1.15 fix bug on ant type replacement by rinex header
*                            fix bug on combined filter for moving-base mode
*           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
*                            add function to read satellite fcb
*                            add function to read stec and troposphere file
*                            add keyword replacement in dcb, erp and ionos file
*           2015/11/13  1.17 add support of L5 antenna phase center paramters
*                            add *.stec and *.trp file for ppp correction
*           2015/11/26  1.18 support opt->freqopt(disable L2)
*           2016/01/12  1.19 add carrier-phase bias correction by ssr
*           2016/07/31  1.20 fix error message problem in rnx2rtkp
*           2016/08/29  1.21 suppress warnings
*           2016/10/10  1.22 fix bug on identification of file fopt->blq
*           2017/06/13  1.23 add smoother of velocity solution
*           2020/11/30  1.24 use API sat2freq() to get carrier frequency
*                            fix bug on select best solution in static mode
*                            delete function to use L2 instead of L5 PCV
*                            writing solution file in binary mode
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss={0};        /* receiver antenna parameters */                   // 接收机天线相位改正参数
static pcvs_t pcvsr={0};        /* satellite antenna parameters */                  // 卫星天线相位改正参数
static obs_t obss={0};          /* observation data */                              // 观测值
static nav_t navs={0};          /* navigation data */                               // 导航电文
static sbs_t sbss={0};          /* sbas messages */
static sta_t stas[MAXRCV];      /* station infomation */
static int nepoch=0;            /* number of observation epochs */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int isbs  =0;            /* current sbas message index */
static int revs  =0;            /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */
static char rtcm_file[1024]=""; /* rtcm data file */
static char rtcm_path[1024]=""; /* rtcm data path */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE *fp_rtcm=NULL;      /* rtcm data file pointer */

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);

    // 输出基准站流动站名
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}

// 输出指定格式的基准站坐标
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];
    const char *sep=opt->sep;
    
    trace(3,"outrpos :\n");
    
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1,5);
            deg2dms(pos[1]*R2D,dms2,5);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}

// 输出结果文件的文件头
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];
    
    trace(3,"outheader: n=%d\n",n);
    
    // 如果是 NMEA、STAT，不需要输出文件头，直接返回
    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }

    // 输出结果文件开头的程序信息
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss.data[i].time;
        te=obss.data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }

    // 相对定位模式调用 outrpos() 输出基准站坐标
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);
    
    // 调用 outsolhead() 输出各结果字段的含义
    outsolhead(fp,sopt);
}

// 寻找下一个时刻的观测值下标，从 i 开始向后的 n 个观测值
// obs->data 的元素已经用sortobs()，根据 time, rcv, sat 排序、去重了
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    // 遍历 obs ，一直正向 i++，找 rcv 对应的观测值，从 i 开始向前的 n 个观测值
    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;

    // 在 i 的基础上加 n++，遍历 obs ，直到流动站变了或时间差大于 DTTOL，找到到这一个时刻观测值的结尾
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);   // 求 i+n 位数据与i数据的时间差 tt
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;          // 时间不同或 rcv 不同，则结束循环
    }
    return n;
}

// 与 nextobsf() 相反，寻找上一个时刻的观测值下标
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}

// 传入时间更新 rtcm.ssr 
// SRR 状态空间域（State Space Representation,SSR ）增强服务，常用于实时 PPP，需要实时获取
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time)
{
    char path[1024];
    int i;
    
    /* open or swap rtcm file */
    reppath(rtcm_file,path,time,"","");
    
    if (strcmp(path,rtcm_path)) {
        strcpy(rtcm_path,path);
        
        if (fp_rtcm) fclose(fp_rtcm);
        fp_rtcm=fopen(path,"rb");
        if (fp_rtcm) {
            rtcm.time=time;
            input_rtcm3f(&rtcm,fp_rtcm);
            trace(2,"rtcm file open: %s\n",path);
        }
    }
    if (!fp_rtcm) return;
    
    /* read rtcm file until current time */
    while (timediff(rtcm.time,time)<1E-3) {
        if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;
        
        /* update ssr corrections */
        for (i=0;i<MAXSAT;i++) {
            if (!rtcm.ssr[i].update||
                rtcm.ssr[i].iod[0]!=rtcm.ssr[i].iod[1]||
                timediff(time,rtcm.ssr[i].t0[0])<-1E-3) continue;
            navs.ssr[i]=rtcm.ssr[i];
            rtcm.ssr[i].update=0;
        }
    }
}

// 取下一个历元的观测数据到 obss，先流动站再基准站，共 n 个
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    int i,  
        nu,     // 当前历元流动站观测值数
        nr      // 当前历元基准站观测值数
        ,n=0;   // 观测值下标
    
    // iobsu ：流动站当前历元索引
    // iobsr ：基准站当前历元索引
    // isbs  ：SBAS信息索引
    // revs  ：0:forward 1:backward

    trace(3,"infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n",revs,iobsu,iobsr,isbs);
    
    if (0<=iobsu&&iobsu<obss.n) {
        settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }

    // 如果是前向滤波
    if (!revs) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        }
        nr=nextobsf(&obss,&iobsr,2);
        if (nr<=0) {
            nr=nextobsf(&obss,&iobsr,2);
        }
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];     // 循环 nu 次，把流动站同一时间、接收机不同卫星的数据加入obs[]
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];     // 循环 nr 次，把基准站的数据加入obs[]
        iobsu+=nu;
        
        /* update sbas corrections */
        while (isbs<sbss.n) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update rtcm ssr corrections */
        if (*rtcm_file) {
            update_rtcm_ssr(obs[0].time);
        }
    }

    // 如果是反向滤波
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        }
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;
        
        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
    }

    // 返回当前历元基准站、流动站观测值总数
    return n;
}

// ssr 载波相位改正
/* carrier-phase bias correction by ssr --------------------------------------*/
static void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
{
    double freq;
    uint8_t code;
    int i,j;

    // 遍历每个观测值每个频率
    for (i=0;i<n;i++) for (j=0;j<NFREQ;j++) {
        code=obs[i].code[j];
        
        if ((freq=sat2freq(obs[i].sat,code,nav))==0.0) continue;
        
        /* correct phase bias (cyc) */
        obs[i].L[j]-=nav->ssr[obs[i].sat-1].pbias[code-1]*freq/CLIGHT;
    }
}

// 从这个函数开始正式对整个流动站和基准站的观测文件内的数据进行逐历元处理
// 每次循环都通过inputobs () 函数读取一个历元的数据，并调用 rtkpos() 函数对该历元的数据进行解算
// 如果是哪个历元出现解算错误，一般都从这个函数开始设置断点向后进行调试
/* process positioning -------------------------------------------------------
 * args:FILE *fp	   		   I/O 输出结果文件指针
 *      const prcopt_t *popt    I   处理选项结构体
 *      const solopt_t *sopt    I   结果选项结构体
 *      const filopt_t *fopt    I   文件选项结构体
 *      int mode			    I   0：forward/backward、1：combined
 ----------------------------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                    int mode)
{
    gtime_t time={0};
    sol_t sol={{0}};
    rtk_t rtk;
    obsd_t obs[MAXOBS*2]; /* for rover and base */
    double rb[3]={0};
    int i,nobs,n,solstatic,pri[]={6,1,2,3,4,5,1,6};
    
    trace(3,"procpos : mode=%d\n",mode);
    
    // 先判断结果是否为静态，处理选项和结果选项都为静态才算静态
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    // 初始化rtk_t，主要是将处理选项结构体 popt 赋值给 rtk 的部分成员
    rtkinit(&rtk,popt);
    rtcm_path[0]='\0';
    
    // 对每一个历元进行遍历求解和输出
    // 获取当前历元观测值数 nobs 以及当前历元各观测记录 obs[MAXOBS*2]
    while ((nobs=inputobs(obs,rtk.sol.stat,popt))>=0) {
        
        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            // 排除禁用卫星的观测值
            // satsys()：传入连续的卫星编号 satellite number
            // 返回卫星系统(SYS_GPS,SYS_GLO,...) ，通过传入的指针 prn 传出 PRN 码
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }
        if (n<=0) continue;
        
        // 如果 ppp 模式设置了 fractional cycle bias 相位的小数轴偏差
        /* carrier-phase bias correction */
        if (!strstr(popt->pppopt,"-ENA_FCB")) {
            corr_phase_bias_ssr(obs,n,&navs);
        }
        // 调用 rtkpos() 进行解算
        if (!rtkpos(&rtk,obs,n,&navs)) continue;
        
        // 单 forward/backward 模式
        if (mode==0) { /* forward/backward */
            // 不是静态模式就直接输出结果
            if (!solstatic) {
                outsol(fp,&rtk.sol,rtk.rb,sopt);
            }
            else if (time.time==0||pri[rtk.sol.stat]<=pri[sol.stat]) {
                sol=rtk.sol;
                for (i=0;i<3;i++) rb[i]=rtk.rb[i];
                if (time.time==0||timediff(rtk.sol.time,time)<0.0) {
                    time=rtk.sol.time;  // 记录上一历元的时间
                }
            }
        }
        else if (!revs) { /* combined-forward */
            if (isolf>=nepoch) return;
            solf[isolf]=rtk.sol;
            for (i=0;i<3;i++) rbf[i+isolf*3]=rtk.rb[i];
            isolf++;
        }
        else { /* combined-backward */
            if (isolb>=nepoch) return;
            solb[isolb]=rtk.sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk.rb[i];
            isolb++;
        }
    }
    if (mode==0&&solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }
    rtkfree(&rtk);
}

// 判断 combine() 结果的有效性，四倍标准差以内有效
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb)
{
    double dr[3],var[3];
    int i;
    char tstr[32];
    
    trace(3,"valcomb :\n");
    
    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];      // 坐标值差 dr 为两坐标相减
        var[i]=solf->qr[i]+solb->qr[i];     // 方差 car 为两相加
    }
    // dr 在限差四倍标准差之内，就合格 return 1，否则 return 0
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */
        
        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}

// 加权平均合并前后向滤波的结果
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0},rr_f[3],rr_b[3],rr_s[3];
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};
    
    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);
    
    // 判断静态模式，处理选项和结果选项都得为静态
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    // i：从前到后，取前向滤波的结果
    // j：从后到前，取后向滤波的结果
    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        // 判断前后向滤波结果的时间差，时间差大于 DTTOL
        // sols、rbs 取时间早的结果，另一个结果的下标不变，进行下一次循环的判断
        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {      // 如果前向时间迟于后向时间
            sols=solf[i];
            // 把前向基站坐标赋值给 rbs[]
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;    // 每次循环都要 j++，这里 j--，相当于保持不变
        }
        // 如果前向时间早于后向时间
        else if (tt>DTTOL) {
            sols=solb[j];
            // 把后向基站坐标赋值给 rbs[]
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;        // i 不变
        }

        // 时间差很小，solution status 不同，sols、rbs 取 solution status 小的结果
        else if (solf[i].stat<solb[j].stat) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (solf[i].stat>solb[j].stat) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        // 时间差很小，solution status 相同
        else {
            sols=solf[i];                           // sols 取前向滤波结果
            sols.time=timeadd(sols.time,-tt/2.0);   // 时间取前后向时间的平均
            
            // 相对定位模式，若结果为固定解，调用 valcomb() 检验，如果失败将 fix 降级为 float
            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {
                
                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j)) sols.stat=SOLQ_FLOAT;
            }

            // 赋值前后向协方差给 Qf、Qb，k+k*3 是取对角线元素
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];  // 赋值非对角线元素
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];
            
            // 调用 smoother() 进行前后向滤波结果结合，位置存在 sols.rr[]，方差存在 sols.qr[]
            if (popt->mode==PMODE_MOVEB) {
                for (k=0;k<3;k++) rr_f[k]=solf[i].rr[k]-rbf[k+i*3];
                for (k=0;k<3;k++) rr_b[k]=solb[j].rr[k]-rbb[k+j*3];
                if (smoother(rr_f,Qf,rr_b,Qb,3,rr_s,Qs)) continue;
                for (k=0;k<3;k++) sols.rr[k]=rbs[k]+rr_s[k];
            }
            else {
                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            }
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];
            
            /* smoother for velocity solution */
            if (popt->dynamics) {
                for (k=0;k<3;k++) {
                    Qf[k+k*3]=solf[i].qv[k];
                    Qb[k+k*3]=solb[j].qv[k];
                }
                Qf[1]=Qf[3]=solf[i].qv[3];
                Qf[5]=Qf[7]=solf[i].qv[4];
                Qf[2]=Qf[6]=solf[i].qv[5];
                Qb[1]=Qb[3]=solb[j].qv[3];
                Qb[5]=Qb[7]=solb[j].qv[4];
                Qb[2]=Qb[6]=solb[j].qv[5];
                if (smoother(solf[i].rr+3,Qf,solb[j].rr+3,Qb,3,sols.rr+3,Qs)) continue;
                sols.qv[0]=(float)Qs[0];
                sols.qv[1]=(float)Qs[4];
                sols.qv[2]=(float)Qs[8];
                sols.qv[3]=(float)Qs[1];
                sols.qv[4]=(float)Qs[5];
                sols.qv[5]=(float)Qs[2];
            }
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,sopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }

    // 循环处理完之后，如果是静态模式且时间存在，调用 outsol() 输出结果
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }
}

// 读取精密星历、SBAS 改正、TEC 格网、初始化 RTCM 解析
/* read prec ephemeris, sbas data, tec grid and open rtcm --------------------*/
static void readpreceph(char **infile, int n, const prcopt_t *prcopt,
                        nav_t *nav, sbs_t *sbs)
{
    seph_t seph0={0};
    int i;
    char *ext;
    
    trace(2,"readpreceph: n=%d\n",n);
    
    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    sbs->n =sbs->nmax =0;
    
    /* read precise ephemeris files */      // 调用 readsp3() 读精密星历 sp3
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */          // 调用 readrnxc() 读精密钟差
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readrnxc(infile[i],nav);
    }
    /* read sbas message files */           // 调用 sbsreadmsg() 读取 SBAS 改正
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        sbsreadmsg(infile[i],prcopt->sbassatsel,sbs);
    }
    /* allocate sbas ephemeris */           // 给 SBAS 星历开辟空间
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;
    
    /* set rtcm file and initialize rtcm struct */  
    rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;
    
    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
            init_rtcm(&rtcm);
            break;
        }
    }
}


/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs)
{
    int i;
    
    trace(3,"freepreceph:\n");
    
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;
    
    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
}

/* read obs and nav data -----------------------------------------------------
 * args:gtime_t ts                  解算开始时间
 *      gtime_t te                  解算结束时间
 *      double ti                   解算时间间隔
 *      char **infile               传入文件路径数组
 *      const int *index            对应文件下标
 *      int n                       infile[] 元素个数
 *      const prcopt_t *prcopt      处理选项
 *      obs_t *obs                  存观测数据 OBS
 *      nav_t *nav                  存导航电文数据 NAV
 *      sta_t *sta                  测站结构体，存观测文件头读取到的一部分的信息
 -----------------------------------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;
    
    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);
    
    // 初始化对所有数据指针
    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;
    
    // 遍历 infile[]，调用 readrnxt() 读取文件
    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;
        
        if (index[i]!=ind) {                // 如果下标和上一次循环的不同
            if (obs->n>nobs) rcv++;         // rcv=1:rover,2:reference
            ind=index[i];                   // 记录当前 index[i] 值到 ind
            nobs=obs->n;                    // 记录观测数据数到 obs->n
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"\n");
        return 0;
    }

    // 调用 sortobs()，根据 time, rcv, sat 
    // 对 obs->data 的元素进行排序、去重，并且得到历元数 nepoch
    /* sort observation data */
    nepoch=sortobs(obs);
    
    // 调用 uniqnav()，进行星历数据的排序去重
    /* delete duplicated ephemeris */
    uniqnav(nav);
    
    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}

// 释放 obs 和 nav
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");
    
    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}

// 通过 nav 和多个 obs 单点定位计算位置，存到 ra[] 中
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];
    
    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);
    
    for (i=0;i<3;i++) ra[i]=0.0;
    
    // for 循环调用 nextobsf() 每次取一个历元观测数据
    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {
        
        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */
        
        // 单点定位，结果存到 sol，再加到 ra[]
        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) continue;
        
        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }

    // ra 除以历元数，得到平均位置
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}

// 从文件中读取测站坐标
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, char *name, double *r)
{
    FILE *fp;
    char buff[256],sname[256],*p,*q;
    double pos[3];
    
    trace(3,"getstapos: file=%s name=%s\n",file,name);
    
    // 以读的方式打开 file
    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    // 循环读取，每次读一行数据，到 \n 或者 256 位结束
    while (fgets(buff,sizeof(buff),fp)) {
        // 如果在行中找到%，截断，赋值\0
        if ((p=strchr(buff,'%'))) *p='\0';
        // 格式化读取，测站位置存到 pos[3]，测站名存到 sname
        if (sscanf(buff,"%lf %lf %lf %s",pos,pos+1,pos+2,sname)<4) continue;
        // 逐字符转大写比较 name、sname
        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            pos[0]*=D2R;
            pos[1]*=D2R;
            pos2ecef(pos,r);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}

// 天线相位中心位置，基准站位置
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;
    
    trace(3,"antpos  : rcvno=%d\n",rcvno);
    
    // 利用基准站的观测文件计算其SPP定位结果作为基准站的坐标
    if (postype==POSOPT_SINGLE) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    // 从 pos 文件读取基准站坐标
    else if (postype==POSOPT_FILE) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
        if (!getstapos(posfile,name,rr)) {
            showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
    }
    // 从基准站的 OBS 观测文件的文件头部分读取基准站坐标
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        // 如果没有坐标数据，报错
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {          
            showmsg("error : no position in rinex header");
            trace(1,"no position position in rinex header\n");
            return 0;
        }
        // 天线相位中心偏差改正
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}

// 
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"openses :\n");
    
    // 调用 readpcv() 读取文件天线相位改正文件
    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
        showmsg("error : no sat ant pcv in %s",fopt->satantp);
        trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    // 调用 readpcv() 读取接收机天线相位改正文件
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
        showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    // 调用 opengeoid() 读取大地水准面数据
    /* open geoid data */
    if (sopt->geoid>0&&*fopt->geoid) {
        if (!opengeoid(sopt->geoid,fopt->geoid)) {
            showmsg("error : no geoid data %s",fopt->geoid);
            trace(2,"no geoid data %s\n",fopt->geoid);
        }
    }
    return 1;
}
// 
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");
    
    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;
    
    /* close geoid data */
    closegeoid();
    
    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;
    
    /* close solution statistics and debug trace */
    rtkclosestat();
    traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv,pcv0={0};
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];
    
    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        nav->pcvs[i]=pcv0;
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(3,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        popt->pcvr[i]=pcv0;
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    
    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}

// 创建结果文件，输出结果文件的文件头
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt)
{
    FILE *fp=stdout;
    
    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);
    
    if (*outfile) {
        createdir(outfile);
        
        if (!(fp=fopen(outfile,"wb"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
    outheader(fp,infile,n,popt,sopt);
    
    if (*outfile) fclose(fp);
    
    return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    
    // ab：以追加的方式打开二进制文件
    return !*outfile?stdout:fopen(outfile,"ab");    
}

//读取各种文件，并将文件中的内容赋值到程序的结构体内。然后，最重要的就是计算基准站的位置以及选择解算类型。
/* execute processing session ------------------------------------------------
 * args:gtime_t ts              I   处理的起始时间，写0表示不限制
 *      gtime_t te              I   处理的起始时间，写0表示不限制
 *      double ti               I   处理的间隔时间 (s)，写0表示不限制，全处理
 *      const prcopt_t *popt    I   处理选项结构体
 *      const solopt_t *sopt    I   结果选项结构体
 *      const filopt_t *fopt    I   文件选项结构体
 *      int flag                I   用于控制输出
 *      char **infile           I   传入文件路径数组首地址
 *      const int *index        I   传入文件路径数组首地址
 *      int n                   I   传入文件数量
 *      char *outfile           I   输出文件的路径，写0表示stdout终端
 * ---------------------------------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                   const solopt_t *sopt, const filopt_t *fopt, int flag,
                   char **infile, const int *index, int n, char *outfile)
{
    FILE *fp;
    prcopt_t popt_=*popt;
    char tracefile[1024],statfile[1024],path[1024],*ext;
    
    trace(3,"execses : n=%d outfile=%s\n",n,outfile);
    
    // 创建 trace 文件，并设置 trace 等级
    /* open debug trace */
    if (flag&&sopt->trace>0) {
        if (*outfile) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
        }
        else {
            strcpy(tracefile,fopt->trace);
        }
        traceclose();
        traceopen(tracefile);
        tracelevel(sopt->trace);
    }
    // 调用 readtec() 读取电离层 TEC 文件，TEC:Total electronic content 总电子含量
    /* read ionosphere data file */
    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I')) {
            reppath(fopt->iono,path,ts,"","");
            readtec(path,&navs,1);
        }
    }

    // 读取地球自转参数 ERP 文件
    /* read erp data */
    if (*fopt->eop) {
        free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
        reppath(fopt->eop,path,ts,"","");
        if (!readerp(path,&navs.erp)) {
            showmsg("error : no erp data %s",path);
            trace(2,"no erp data %s\n",path);
        }
    }

    // 调用 readobsnav() 读取 OBS 和 NAV 文件
    /* read obs and nav data */
    if (!readobsnav(ts,te,ti,infile,index,n,&popt_,&obss,&navs,stas)) return 0;
    
    // 调用 readdcb() 读取差分码偏差 DCB 参数
    /* read dcb parameters */
    if (*fopt->dcb) {
        reppath(fopt->dcb,path,ts,"","");
        readdcb(path,&navs,stas);
    }

    // 调用 setpcv() 设置天线相位中心变化参数 PCV 
    /* set antenna paramters */
    if (popt_.mode!=PMODE_SINGLE) {
        setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
               stas);
    }

    // 调用 readotl() 读取潮汐参数
    /* read ocean tide loading parameters */
    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
        readotl(&popt_,fopt->blq,stas);
    }
    // FIXED 模式，调用 antpos() 得到流动站坐标
    /* rover/reference fixed position */
    if (popt_.mode==PMODE_FIXED) {
        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
    }
    // DGPS、KINEMA、STATIC 模式，调用 antpos() 得到基准站坐标
    else if (PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC) {
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
    }

    // 创建解算状态文件（内容包括残差、高度角等，解算效果不好的时候输出来看具体为啥算不好）
    /* open solution statistics */
    if (flag&&sopt->sstat>0) {
        strcpy(statfile,outfile);
        strcat(statfile,".stat");
        rtkclosestat();
        rtkopenstat(statfile,sopt->sstat);
    }

    // 调用 outhead() 创建结果文件，输出结果文件的文件头
    /* write header to output file */
    if (flag&&!outhead(outfile,infile,n,&popt_,sopt)) {
        freeobsnav(&obss,&navs);
        return 0;
    }
    iobsu=iobsr=isbs=revs=aborts=0;
    
    // 针对不同的滤波解算类型，inputobs() 函数内读取文件数据的顺序不同
    if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) {
        if ((fp=openfile(outfile))) {
            procpos(fp,&popt_,sopt,0); /* forward */    // 前向滤波
            fclose(fp);
        }
    }
    else if (popt_.soltype==1) {
        if ((fp=openfile(outfile))) {
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;
            procpos(fp,&popt_,sopt,0); /* backward */   // 后向滤波
            fclose(fp);
        }
    }
    else { /* combined */
        // 开辟内存空间
        solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);     // 前向结果
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);     // 后向结果
        rbf=(double *)malloc(sizeof(double)*nepoch*3);  // 前向基准站坐标
        rbb=(double *)malloc(sizeof(double)*nepoch*3);  // 后向基准站坐标
        
        if (solf&&solb) {   // 判断内存开辟成功
            isolf=isolb=0;
            procpos(NULL,&popt_,sopt,1); /* forward */      // 前向滤波
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;    
            procpos(NULL,&popt_,sopt,1); /* backward */     // 后向滤波
            
            /* combine forward/backward solutions */
            if (!aborts&&(fp=openfile(outfile))) {
                combres(fp,&popt_,sopt);
                fclose(fp);
            }
        }
        else showmsg("error : memory allocation");
        free(solf);
        free(solb);
        free(rbf);
        free(rbb);
    }
    /* free obs and nav data */
    freeobsnav(&obss,&navs);
    
    return aborts?1:0;
}

//对每个流动站执行定位过程
/* execute processing session for each rover ---------------------------------
 * args:gtime_t ts              I   处理的起始时间，写0表示不限制
 *      gtime_t te              I   处理的起始时间，写0表示不限制
 *      double ti               I   处理的间隔时间 (s)，写0表示不限制，全处理
 *      const prcopt_t *popt    I   处理选项结构体
 *      const solopt_t *sopt    I   结果选项结构体
 *      const filopt_t *fopt    I   文件选项结构体
 *      int flag                I   用于控制输出
 *      char **infile           I   传入文件路径数组首地址
 *      const int *index        I   传入文件路径数组首地址
 *      int n                   I   传入文件数量
 *      char *outfile           I   输出文件的路径，写0表示stdout终端
 *      const char *rov         I   流动站ID列表，空格隔开
 -----------------------------------------------------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*rov_,*p,*q,s[64]="";
    
    trace(3,"execses_r: n=%d outfile=%s\n",n,outfile);
    
    for (i=0;i<n;i++) if (strstr(infile[i],"%r")) break;
    
    // 如果某个 infile[i] 含有流动站 ID 的替换符
    if (i<n) { /* include rover keywords */
        if (!(rov_=(char *)malloc(strlen(rov)+1))) return 0;
        strcpy(rov_,rov);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(rov_); for (;i>=0;i--) free(ifile[i]);
                return 0;
            }
        }
        for (p=rov_;;p=q+1) { /* for each rover */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_rov,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,p,"");
                reppath(outfile,ofile,t0,p,"");
                
                /* execute processing session */
                stat=execses(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile);
            }
            if (stat==1||!q) break;
        }
        free(rov_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        /* execute processing session */
        stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    }
    return stat;
}

//对每个基准站执行定位过程
/* execute processing session for each base station --------------------------
 * args:gtime_t ts              I   处理的起始时间，写0表示不限制
 *      gtime_t te              I   处理的起始时间，写0表示不限制
 *      double ti               I   处理的间隔时间 (s)，写0表示不限制，全处理
 *      const prcopt_t *popt    I   处理选项结构体
 *      const solopt_t *sopt    I   结果选项结构体
 *      const filopt_t *fopt    I   文件选项结构体
 *      int flag                I   用于控制输出
 *      char **infile           I   传入文件路径数组首地址
 *      const int *index        I   传入文件路径数组首地址
 *      int n                   I   传入文件数量
 *      char *outfile           I   输出文件的路径，写 0 表示输出到 stdout 终端
 *      const char *rov         I   流动站ID列表，空格隔开
 *      const char *base        I   基准站ID列表，空格隔开
 -----------------------------------------------------------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov, const char *base)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*base_,*p,*q,s[64];
    
    trace(3,"execses_b: n=%d outfile=%s\n",n,outfile);
    
    /* read prec ephemeris and sbas data */
    // 调用 readpreceph() 读取精密星历和 SBAS 数据
    readpreceph(infile,n,popt,&navs,&sbss);
    
    // %b：基准站 ID 的替换符
    for (i=0;i<n;i++) if (strstr(infile[i],"%b")) break;
    
    // 如果某个 infile[i] 含有基准站 ID 的替换符，找基准站文件
    if (i<n) { /* include base station keywords */
        // 为 base_ 开辟空间，将 base 赋值给 base_
        if (!(base_=(char *)malloc(strlen(base)+1))) {
            freepreceph(&navs,&sbss);
            return 0;
        }
        strcpy(base_,base);
        
        // 为 ifile[] 开辟空间
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(base_); for (;i>=0;i--) free(ifile[i]);
                freepreceph(&navs,&sbss);
                return 0;
            }
        }

        // 遍历 base_ 基准站字符串
        for (p=base_;;p=q+1) { /* for each base station */
            if ((q=strchr(p,' '))) *q='\0';     // 拆出一个基准站
            
            if (*p) {
                strcpy(proc_base,p);            // 把基准站名赋值给 proc_base
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                // 循环替换 infile[i] 里的基准站 ID 的替换符到 ifile[i]
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,"",p);
                // 替换 outfile 里的基准站 ID 替换符到 ofile
                reppath(outfile,ofile,t0,"",p);
                // 调用 execses_r() 进行下一步解算
                stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile,rov);
            }
            if (stat==1||!q) break;
        }
        free(base_); for (i=0;i<n;i++) free(ifile[i]);
    }
    // infile[i] 都没有有基准站 ID 的替换符，直接调用 execses_r() 进行下一步解算
    else {
        stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);
    }
    /* free prec ephemeris and sbas data */
    freepreceph(&navs,&sbss);
    
    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise 
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first 
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile,
                   const char *rov, const char *base)
{
    gtime_t tts,                // 解算开始时间
            tte,                // 解算结束时间
            ttte;               // 读取星历文件的结束时间
    double tunit,               // 
           tss;                 // 
    int i,j,k,                  // 循环和数组下标控制
        nf,                     // 文件路径数组下标控制
        stat=0,                 // 接收返回状态值，为1
        week,                   // 用于存 GPST 的周
        flag=1,                 // 
        index[MAXINFILE]={0};   // 
    char *ifile[MAXINFILE],     // 
         ofile[1024],           // 
         *ext;                  // 
    
    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);
    
    // 调用 openses() 开始处理，文件读取，赋值 navs、pcvs、pcvsr
    /* open processing session */
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;
    
    // 判断 ts、te、tint 时间是否大于 0，即是否设置了合理的时间，有三种情况
    // 1. if (ts.time!=0&&te.time!=0&&tu>=0.0) {...
    // 2. else if (ts.time!=0) {...
    // 3. else...                                  
    
    if (ts.time!=0&&te.time!=0&&tu>=0.0) {
        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
            closeses(&navs,&pcvss,&pcvsr);
            return 0;
        }
        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }
        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
        settspan(ts,te);
        tunit=tu<86400.0?tu:86400.0;                        // tunit：如果 tu 小于一天就为 tu，否则为一天
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);
        
        // 根据解算时间单元，分时间段循环处理，算出来 tts>te 或过程有错误，结束循环
        // 很多时候解算单元时间直接设 0.0，只循环一次，tts=ts，tte=te
        for (i=0;;i++) { /* for each periods */
            tts=gpst2time(week,tss+i*tu);       // 解算单元开始时间，每次循环加上一个 i 个 tu
            tte=timeadd(tts,tu-DTTOL);          // 解算结束时间 tte=tu-DTTOL
            if (timediff(tts,te)>0.0) break;    // 算出来 tts>te 结束循环
            if (timediff(tts,ts)<0.0) tts=ts;   // 分时间段后 tts 若早于 ts，设为 ts
            if (timediff(tte,te)>0.0) tte=te;   // 分时间段后 tte 若早于 te，设为 te
            
            strcpy(proc_rov ,"");
            strcpy(proc_base,"");
            if (checkbrk("reading    : %s",time_str(tts,0))) {
                stat=1;
                break;
            }
            for (j=k=nf=0;j<n;j++) {
                
                ext=strrchr(infile[j],'.');
                
                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
                    strcpy(ifile[nf++],infile[j]);
                }
                else {
                    // 星历文件，包括精密星历和广播星历
                    /* include next day precise ephemeris or rinex brdc nav */
                    ttte=tte;
                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                        ttte=timeadd(ttte,3600.0);          // 精密星历加一小时
                    }
                    else if (strstr(infile[j],"brdc")) {
                        ttte=timeadd(ttte,7200.0);          // 广播星历加两小时
                    }
                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
                }
                while (k<nf) index[k++]=j;
                
                if (nf>=MAXINFILE) {
                    trace(2,"too many input files. trancated\n");
                    break;
                }
            }
            if (!reppath(outfile,ofile,tts,"","")&&i>0) flag=0;
            
            /* execute processing session */
            stat=execses_b(tts,tte,ti,popt,sopt,fopt,flag,ifile,index,nf,ofile,
                           rov,base);
            
            if (stat==1) break;
        }
        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }
    else if (ts.time!=0) {
        for (i=0;i<n&&i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                return -1;
            }
            reppath(infile[i],ifile[i],ts,"","");
            index[i]=i;
        }
        reppath(outfile,ofile,ts,"","");
        
        /* execute processing session */
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,ifile,index,n,ofile,rov,
                       base);
        
        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
    }
    else {
        for (i=0;i<n;i++) index[i]=i;
        
        /* execute processing session */
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,
                       base);
    }
    /* close processing session */
    closeses(&navs,&pcvss,&pcvsr);
    
    return stat;
}
