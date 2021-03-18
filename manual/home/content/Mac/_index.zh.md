---
title: "Mac"
date: 2018-12-29T11:02:05+06:00
icon: "ti-pie-chart"
logo: "images/mac.svg"
appSnapShot: "../images/sweos_mac_zh.png"
systemVersion: Catalina
installerFormats: ["zip", "dmg"]
warning: 由于没有经费支持苹果开发者账户($99/年)，故swEOS应用程序没有被Apple开发中心公证，这会导致在安装使用的时候出现一些问题，解决方案如下所示！
description: "桌面应用程序、纯终端命令程序、C++动态库"
type : "apps"
---

#### 如何安装没有经过Apple公证的软件

#### 1. 使用没有被公证的软件会出现什么状况 ?

macOS从Catalina版本以后对安全性有了进一步提升，对软件的检查更严格。
系统中引入了看门狗(Gatekeeper😄)技术(更多详情见[Apple support](https://support.apple.com/en-us/HT202491))，用户在第一次运行软件的时候看门狗就会检查这个软件是否被苹果的软件中心公证过。
如果没有被公证，则提示用户在安全中心中手动允许运行。
一般情况下，就在系统的`安全与隐私`窗口中点击允许运行即可，但是`swEOS`的可视化系统是基于[Qt](https://www.qt.io)和[VTK](https://vtk.org)，并行计算时基于[OpenMP](https://www.openmp.org)的，所以有几十个依赖库都需要手动允许运行。
一个一个点按钮，显然是要累死人的节奏啊，不过好消息是我已经找到了解决方案，见下面的两个步骤。

<image src=../images/app_mac_install4.png width=100% ></image>
<image src=../images/app_mac_install2.jpg width=100% ></image>

#### 2. 解决方案

**2.1 第一步** : 将`swEOS`拖进`Applications`文件夹 (见下图). 如果刚才你已经完成了这一步，则跳过此步骤。

<image src=../images/app_mac_install3.png width=100% ></image>

**2.2 Second step** : 在终端运行命令通知看门狗swEOS是可信的

有两种方式让看门狗听话，允许经过主人认可的软件进门😎，很容易，只需要用管理员权限(`sudo`)在终端下运行一行命令即可。

* 第一种方式: 告诉看门狗只允许`swEOS.app`，这种方式是我推荐的，更安全一些。

```python
sudo xattr -r -d com.apple.quarantine /Applications/swEOS.app
```

* 第二种方式: 告诉看门狗允许所有的app都可以运行，这显然是存在一定安全隐患的，这里只做介绍但并不推荐使用。

```python
sudo spctl --master-disable
```