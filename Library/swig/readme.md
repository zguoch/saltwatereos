
## Build API

### Dependence

- swig
- Python
- TCL
- cmake 

### Mac and Linux
There are switches of Python, TCL in the main [CMakeLists.txt](../CMakeLists.txt). For Mac and Linux, it is pretty easy. 

### Windows
The compiler can be Microsoft Visual Studio (e.g. 2017 community). The python version can be anaconda 3. swig and TCL can be installed from official cite binary package.
Because the environment variable or cmake package variable is not clear on windows system, just use the CMake GUI version to configure the project, the related path of swig, python and TCL can be set following the GUI steps.
**Note** that if the python is 64bit, the project has to be configured as **x64** and **release** type.

# JavaScript 

The latest version of nodejs (>v8) seems doesn't support swig, at least can not correctly build. Therefore, we have to use the lower version of node. The version 8 works!
A pretty nice way to manage versions of node is the tool of `nvm`, it can be install from command of `brew install nvm` on mac, and then do a little bit configuration,

```
export NVM_DIR="$HOME/.nvm"
source $(brew --prefix nvm)/nvm.sh
```
## Using nvm install node v8
Firstly, you have to using command of `nvm ls` to list all the available versions of nodejs, the results looks like this,
```
nvm ls
        v6.14.4
->      v8.17.0
       v14.15.1
        v15.2.1
default -> 8.17.0 (-> v8.17.0)
node -> stable (-> v15.2.1) (default)
stable -> 15.2 (-> v15.2.1) (default)
iojs -> N/A (default)
unstable -> N/A (default)
lts/* -> lts/fermium (-> v14.15.1)
lts/argon -> v4.9.1 (-> N/A)
lts/boron -> v6.17.1 (-> N/A)
lts/carbon -> v8.17.0
lts/dubnium -> v10.23.0 (-> N/A)
lts/erbium -> v12.19.1 (-> N/A)
lts/fermium -> v14.15.1
```

Then install `v8.17.0` using command of `nvm install 8.17.0`, and finally set it to default `nvm alias default 8.17.0`. Now it should works on Mac, check node and npm version,

```
➜  swig git:(master) ✗ node --version
v8.17.0
➜  swig git:(master) ✗ npm --version          
6.13.4
➜  swig git:(master) ✗ which npm
/Users/zguo/.nvm/versions/node/v8.17.0/bin/npm
➜  swig git:(master) ✗ 
```

Using v8 node install `node-gyp`: `npm install node-gyp -g` and then build js API, `node-gyp configure build`