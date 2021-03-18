---
title: "Mac"
date: 2018-12-29T11:02:05+06:00
icon: "ti-panel"
logo: "images/mac.svg"
appSnapShot: "../images/sweos_mac_en.png"
systemVersion: Catalina
installerFormats: ["zip", "dmg"]
warning: The app has not been notarized by Apple because there is no funding to support a Apple Developer ID. This could cause some installation issues, please find the installation tips below.
description: "Desktop app, command line app, c++ library"
type : "apps"
---

#### How to install without notarization

#### 1. What happens when starting to use a non-notarized app ?

macOS includes a technology called Gatekeeper, that's designed to ensure that only trusted software runs on your Mac (see [Apple support](https://support.apple.com/en-us/HT202491) for details). 
That means the non-notarized app will be blocked by the Gatekeeper.
Usually just click `Allow Anyway` then you can use the app(see the second snapshot below), but `swEOS` depends on a lot of third-party libraries, include [Qt](https://www.qt.io), [VTK](https://vtk.org) and [OpenMP](https://www.openmp.org), therefore it is not a easy thing to click dozens of 
`Allow Anyway` button. 
**There is a one-time process can help you to solve this problem, see the following two steps.**

<image src=../images/app_mac_install4.png width=100% ></image>
<image src=../images/app_mac_install2.jpg width=100% ></image>

#### 2. How to solve this issue ?

There are two ways to get rid of this annoying security check information, but I recommend only the first one which only "open the door" for swEOS app. The solution is pretty easy, just run one line of command with **super user** (`sudo`) rights.

**2.1 First step** : drag `swEOS` app to the `Applications` folder (see snapshot below). If you already finished this step before, just skip this step.

<image src=../images/app_mac_install3.png width=100% ></image>

**2.2 Second step** : run command in the terminal

* The first option: tell the Gatekeeper only allow `swEOS.app`, this is a much more secure way.

```python
sudo xattr -r -d com.apple.quarantine /Applications/swEOS.app
```

* The second option: tell the Gatekeeper allow all apps. Obviously there are some security risks in this way. Here I just introduce this option, but dont's recommend to use is actually.

```python
sudo spctl --master-disable
```