rm *.o *.so 
rm *_wrap* example.py

python_inc=/Users/zguo/.pyenv/versions/anaconda3-2019.10/include/python3.7m
# 
python_lib=/Users/zguo/.pyenv/versions/anaconda3-2019.10/lib

# # python 
# swig -c++ -python example.i
# g++ -c example_wrap.cxx example.cpp  -I$python_inc
# # gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
# g++ -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace

# # tcl 
# swig -c++ -tcl example.i
# g++ -c example_wrap.cxx example.cpp
# # gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
# g++ -shared -o example.so  example_wrap.o example.o -ltcl -flat_namespace

# # js 
# # -jsc, -v8 -DV8_VERSION=0x032530
# # swig -c++ -javascript -v8 example.i
# # swig -c++ -javascript -jsc example.i
# # swig -javascript -node -c++ -v8 -DV8_VERSION=0x041027 example.i
# swig -c++ -javascript -node example.i
# g++ -c example_wrap.cxx example.cpp -I/usr/local/Cellar/node/14.4.0//include/node
# # node-gyp configure # build
# # g++ -c example_wrap.cxx example.cpp
# # # gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
# # g++ -shared -o example.so  example_wrap.o example.o -ltcl -flat_namespace

# R
swig -c++ -r example.i
swig -c++ -r -o example_wrap.cpp example.i
# R CMD SHLIB example_wrap.cpp example.cpp
R CMD SHLIB -o example.so example_wrap.cpp example.cpp