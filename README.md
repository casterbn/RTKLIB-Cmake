# RTKLIB-CMake

[TOC]

## 一、项目介绍

### 1、初衷

基于 CMake 构建的 RTKLIB 项目，同时构建 rnx2rtkp、rtkrcv、str2str、convbin、pos2kml 五大命令行程序，并且加上详细的中文注释，方便进行使用和学习，在 Windows 和 Linux 下都可以运行，用 VS、VScode、Clion 都可调试。



![RTKLIB](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/RTKLIB.png)

### 2、相关链接

* **RTKLIB 官网**：https://www.rtklib.com/

* **RTKLIB-2.4.3 源码**：https://github.com/tomojitakasu/RTKLIB/tree/rtklib_2.4.3
* **RTKLIB-2.4.3 程序**：https://github.com/tomojitakasu/RTKLIB_bin/tree/rtklib_2.4.3
* **RTKLIB-demo5 源码**：https://github.com/rtklibexplorer/RTKLIB



### 3、使用说明







### 4、开源许可







## 二、VScode + Linux 下编译调试

### 1、WSL 环境配置

* **安装 WSL**：建议安装 WSL2，先在系统设置里启用虚拟机，然后在微软商城安装 Ubuntu，



* **更新软件列表**

  ```bash
  sudo apt update
  ```

* **安装编译器和调试器** GCC、G++、GDB

  ```bash
  sudo apt install build-essential gdb
  ```

* **安装 CMake**

  ```bash
  sudo apt install cmake
  ```

* **检测安装是否成功**

  ```bash
  cmake -version
  ```

* **安装 ssh**

  ```bash
  sudo apt install ssh
  ```

* 其它推荐软件的安装：安装方式都是 `sudo apt install 软件名`

  - **vim**：LInux 系统必备的文件编辑器，建议没接触过的朋友花半小时学一下操作。
  - **git**：
  - **ranger**：命令行下的资源管理器。
  - **mc**：也是命令行下的资源管理器。
  - **exa**：代替 ls。
  - **tree**：用于查看树状目录结构。
  - **cloc**：统计代码行数。
  - **ack**：替代 grep 查找。
  - **glances**：系统监控工具。
  - **colordiff**：替代 bat。
  - **bat**：替代 cat 可以简单查看文件内容。
  - **dstat**：

### 2、VScode 插件安装

插件下载，可能会比较慢，可以官网下载，然后导入VScode，

注意插件要装在 WSL 中，而不是本地(Local)。

必须要安装的插件包括 C++、CMake、SSH，你也可以装一下别的辅助编程的插件。



### 3、VScode 连接 Linux







### 4、把 RTKLIB 编译成第三方库



* 指定最小 CMake 版本、子项目名

  ```cmake
  cmake_minimum_required(VERSION 3.0)
  project(rtklib)
  ```

* 设置编译时 gcc 参数：

  ```cmake
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3 -ansi -pedantic")
  set(CMAKE_C_FLAGS "-Wno-unused-but-set-variable -Wno-format-overflow -Wno-unused-result -Wpointer-to-int-cast")
  ```

* 指定头文件目录：

  ```cmake
  include_directories(include)
  ```

* 指定可执行文件的输出路径为：rtklib/bin、库文件的输出路径为：rtklib/lib：

  ```cmake
  set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
  ```

* 将 src、src/rcv 目录下源文件加到 DIR_SRCS 列表：

  ```cmake
  aux_source_directory(src DIR_SRCS_RTKLIB)
  aux_source_directory(src/rcv DIR_SRCS_RTKLIB_RCV)
  list(APPEND DIR_SRCS ${DIR_SRCS_RTKLIB} ${DIR_SRCS_RTKLIB_RCV})
  ```

* 把代码编译成动态库，链接上 pthread m 库：

  ```cmake
  add_library(${PROJECT_NAME} SHARED ${DIR_SRCS})
  target_link_libraries(${PROJECT_NAME} pthread m)
  target_include_directories(${PROJECT_NAME}
      PUBLIC ${PROJECT_SOURCE_DIR}/include
  )
  ```

* 如果是 WIN32 还有链接上 wsock32 ws2_32 winmm 库，加上宏定义 -DWIN_DLL：

  ```cmake
  if(WIN32)
    target_link_libraries(${PROJECT_NAME} wsock32 ws2_32 winmm)
    add_definitions(-DWIN_DLL)
  endif()
  ```





### 5、编译调试命令行程序、链接 RTKLIB

RTKLIB APP 目录下有 5 个命令行程序

* **rnx2rtkp**：后处理定位解算
* **rtkrcv**：实时定位解算
* **str2str**：数据流转换播发
* **convbin**：数据转换
* **pos2kml**：定位结果转谷歌地图数据格式





## 三、其它环境的编译调试

> 输命令行参数的几种方式：
>
> * **VS**：
> * **VScode**：
> * **Clion**：
> * **Windows\Linux  终端**：`可执行文件 参数1 参数2 ......`

### 1、Windows + VScode





### 2、CLion + Linux





### 3、CLion + Windows





### 4、VS + Windows





### 5、VS + Linux







---

## 四、rnx2rtkp：基于 RINEX 后处理定位解算

```cmake
rnx2rtkp [option ...] file file [...] 
```

### 1、简介

rnx2rtkp 全称 RINEX to RTK pos，通过原始 RINEX 文件，输出 RTKLIB 的定位坐标，下图可以很好的表示整个程序的逻辑：

![image-20231016123911526](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/image-20231016123911526.png)

* 使用方式：`rnx2rtkp [option]... file file [...] `
* 读取 RINEX：OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS 等文件，计算接收机、流动站坐标，并输出结果。
* 对于相对定位，第一个 OBS 观测值文件需含接收机、流动站观测值，第二个 OBS 文件需含基准站观测值。
* 输入文件至少要有一个星历文件，RINEX NAV/GNAV/HNAV 。
* 想用 SP3 精密星历文件，需提供 .sp3/.eph 文件的路径。
* 输入文件路径可包含通配符 *，为了防止与命令行命令冲突，要用 `"..."`  括起带通配符符路径。



### 2、配置选项

* `filopt_t`：**文件选项**，存结果输出、Trace、各种改正文件路径，不包括星历文件和观测文件。
* `solopt_t`：**结果选项**，可以设置结果输出形式（ENU、ECEF、NMEA、BLH），小数位数，是否输出文件头，是否输出速度等。
* `prcopt_t`：**处理选项**，是配置的重头戏，可以先看 [postpos 的用法](https://www.bilibili.com/video/BV1m5411Y7xV) 学习界面程序的配置方式。写代码配置和界面程序需要配置的东西是一样的，只是从在界面上选，换成在了代码里给对应字段赋值，或者在配置文件中设置。

![](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/RTKLIB%25E9%2585%258D%25E7%25BD%25AE%25E9%2580%2589%25E9%25A1%25B9.png)

### 3、配置文件





### 4、命令行参数

* **-？**：打印 help
* **-k** file：配置文件的输入选项，默认值是 [off]
* **-o** file：输出文件选项，默认值是 [stdout]
* **-ts** ds ts：设置开始解算时间`(ds=y/m/d ts=h:m:s) `，默认值是 [obs start time] 
* **-te** de ds：设置结束解算时间`(de=y/m/d te=h:m:s) `，默认值是 [obs end time] 
* **-ti** tint：设置解算时间间隔频率`(sec) `，默认值是[all]
* **-p** mode：设置解算模式，(**0**:single,**1**:dgps,**2**:kinematic,**3**:static,**4**:moving-base,**5**:fixed,**6**:ppp-kinematic,**7**:ppp-static)，默认值是 [2]
* **-m** mask：设置截止高度角，`(deg) `,默认值是 [15]
* **-sys** s：设置用于计算的导航系统，`(s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) `，默认值是 [G|R] ，想用除 GPS 以外的系统，还得加宏定义 ENAGLO、ENACMP、ENAGAL
* **-f** freq：设置用于计算的频率，` (1:L1,2:L1+L2,3:L1+L2+L5) `，默认值是 [2]
* **-v** thres：设置整周模糊度 Ratio 值，写 0.0 为不固定整周模糊度，默认值是 [3.0] 
* **-b**：后向滤波
* **-c**：前后向滤波组合
* **-i**：单历元模糊度固定 instantaneous 
* **-h**：fix and hold 模糊度固定
* **-e**：输出 XYZ-ecef 坐标
* **-a**：输出 ENU-baseline
* **-n**：输出 NMEA-0183 GGA
* **-g**：输出经纬度格式为 ddd mm ss.ss ，默认为 [ddd.ddd] 
* **-t**：输出时间格式为 yyyy/mm/dd hh:mm:ss.ss ，默认为 [sssss.ss] 
* **-u**：输出为 UTC 时间，默认为 [gpst] 
* **-d** col：设置时间的小数位数，默认为 [3] 
* **-s** sep：设置文件分隔符，要写在单引号中，默认为 [' '] 
* **-r** x y z：基站位置 ECEF-XYZ (m)，默认 [average of single pos] ，流动站位置用于 fixed 模式
* **-l** lat lon hgt：基站位置 LLH (deg/m)，默认 [average of single pos]，流动站位置用于 fixed模式
* **-y** level：输出结果信息 (**0**:off,**1**:states,**2**:residuals) ，默认为 [0] 
* **-x** level：输出 debug trace 等级，默认为 [0] 

### 5、函数调用关系

![image-20231012212121216](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/image-20231012212121216.png)

### 6、执行流程

![image-20231025155540386](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/image-20231025155540386.png)

---

## 五、rtkrcv：实时定位解算

```cmake
rtkrcv [-s][-p port|-d dev][-o file][-t level]
```

### 1、简介

命令行实时定位解算程序。







### 2、命令行参数

*     `-s`：程序启动的时候开启 RTK 定位解算
*     `-p port`：port number for telnet console
*     `-m port`：port number for monitor stream
*     `-d dev`：terminal device for console
*     `-o file`：处理选项文件
*     `-w pwd`：login password for remote console ("": no password)
*     `-r level`：output solution status file (0:off,1:states,2:residuals)
*     `-t level`：debug trace level (0:off,1-5:on)
*     `-sta sta`：station name for receiver dcb

### 3、终端命令

* `start`：开启实时解算，如果程序执行时传入了 -s 参数就不需要再用这个命令。
* `stop`：停止定位解算。
* `restart`：重启定位解算，如果处理选项重新设置了，发送这个命令使能新选项。
* `solution [cycle]`：输出实时定位结果，如果
* `status [cycle]`：输出解算状态，
* `satellite [-n] [cycle]`：输出卫星状态，
* `observ [-n] [cycle]`：
* `navidata [cycle]`：
* `stream [cycle]`：
* `error`：
* `option [opt]`：
* `set opt [val]`：
* `load [file]`：
* `save [file]`：
* `log [file|off]`：
* `help|? [path]`：
* `exit`：
* `shutdown`：
* `!command [arg...]`：

### 4、函数调用关系







### 5、程序执行流程

![image-20231024203128852](https://pic-bed-1316053657.cos.ap-nanjing.myqcloud.com/img/image-20231024203128852.png)

---

## 六、str2str 数据流存储转发

```cmake
str2str -in stream[#...] -out stream[#...] [-out stream[#...]...] [options] 
```

```yaml
OPTIONS 
	-in stream[#format] input stream path and format 
	-out stream[#format] output stream path and format 
```

### 1、简介

从数据流中输入数据，并将其分割和输出到多个数据流中，输入流可以是串行、tcp 客户端、tcp 服务器、ntrip 客户端或文件。输出流可以是串行、tcp 客户端、tcp 服务器、ntrip 服务器或文件。str2str 是常驻应用程序。要停止它：

* 如果运行在前台，则在控制台中键入 ctr-c；
* 如果运行在后台，则向后台进程发送 SIGINT 信号。

如果输入流和输出流都遵循 #format 输入信息的格式将被转换为输出格式。要指定输出使用 -msg 选项。如果省略选项 -in 或 -out，则输入为 stdin，"... "输出为 stdout、输入使用stdin，输出使用 stdout。如果选项 -in 或 -out 中的流为空，也会使用 stdin 或 stdout



### 2、配置文件





### 3、命令行参数

使用方法：`str2str [-in stream] [-out stream [-out stream...]] [options]`



数据流路径 stream path：

```yaml
serial       : serial://port[:brate[:bsize[:parity[:stopb[:fctr]]]]]
tcp server   : tcpsvr://:port
tcp client   : tcpcli://addr[:port]
ntrip client : ntrip://[user[:passwd]@]addr[:port][/mntpnt]
ntrip server : ntrips://[:passwd@]addr[:port]/mntpnt[:str] (only out)
ntrip caster : ntripc://[user:passwd@][:port]/mntpnt[:srctbl] (only out)
file         : [file://]path[::T][::+start][::xseppd][::S=swap]
```

数据格式 format：

```yaml
rtcm2        : RTCM 2 (only in)
rtcm3        : RTCM 3
nov          : NovAtel OEMV/4/6,OEMStar (only in)
oem3         : NovAtel OEM3 (only in)
ubx          : ublox LEA-4T/5T/6T (only in)
ss2          : NovAtel Superstar II (only in)
hemis        : Hemisphere Eclipse/Crescent (only in)
stq          : SkyTraq S1315F (only in)
javad        : Javad (only in)
nvs          : NVS BINR (only in)
binex        : BINEX (only in)
rt17         : Trimble RT17 (only in)
sbf          : Septentrio SBF (only in)
```

选项 option：

```yaml
-sta sta          station id
-opt opt          receiver dependent options
-s  msec          timeout time (ms) [10000]
-r  msec          reconnect interval (ms) [10000]
-n  msec          nmea request cycle (m) [0]
-f  sec           file swap margin (s) [30]
-c  file          input commands file [no]
-c1 file          output 1 commands file [no]
-c2 file          output 2 commands file [no]
-c3 file          output 3 commands file [no]
-c4 file          output 4 commands file [no]
-p  lat lon hgt   station position (latitude/longitude/height) (deg,m)
-px x y z         station position (x/y/z-ecef) (m)
-a  antinfo       antenna info (separated by ,)
-i  rcvinfo       receiver info (separated by ,)
-o  e n u         antenna offset (e,n,u) (m)
-l  local_dir     ftp/http local directory []
-x  proxy_addr    http/ntrip proxy address [no]
-b  str_no        relay back messages from output str to input str [no]
-t  level         trace level [0]
-fl file          log file [str2str.trace]
-h                print help
```







---

## 七、convbin：数据文件格式转换

```yaml
convbin [-ts y/m/d h:m:s] [-te y/m/d h:m:s] [-ti tint] [-r format] [-ro opts] 
 		[-f freq] [-hc comment] [-hm marker] [-hn markno] [-ht marktype] 
	 	[-ho observ] [-hr rec] [-ha ant] [-hp pos] [-hd delta] [-v ver] [-od] 
 		[-os] [-x sat] [-y sys] [-d dir] [-c satid] [-o ofile] [-n nfile] 
 		[-g gfile] [-h hfile] [-q qfile] [-s sfile] file 
```

### 1、简介





### 2、配置文件





### 3、命令行参数





### 4、函数调用关系





### 5、执行流程









---

## 八、pos2kml：定位结果转谷歌地图格式

```cmake
pos2kml [option ...] file [...] 
```

### 1、简介





### 2、命令行参数





### 3、函数调用关系





### 4、执行流程













