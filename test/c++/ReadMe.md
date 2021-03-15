
The swEOS library can be called by any c++ source code, you could try as following.

**Requirement**

1. [CMake](https://cmake.org)
2. C++ compiler

The following test could work in most of macOS and Linux system. For the Windows system with MSVC (e.g. Visual Studio 2017 Community version), just open the generated `.sln` file after running `cmake ..` command, and then build the project.

**Option 1: default settings in this demo**

```bash
mkdir build
cd build
cmake ..
make 
./test_swEOS
```

**Option 2: specify the `SWEOS_DIR` path which contains `include` and `lib` of swEOS library**

```bash
mkdir build
cd build
cmake -DSWEOS_DIR=../../
make
./test_swEOS
```