rm *.o *.so 
rm *_wrap* example.py

python_inc=/Users/zguo/.pyenv/versions/anaconda3-2019.10/include/python3.7m
# 
python_lib=/Users/zguo/.pyenv/versions/anaconda3-2019.10/lib

# python 
swig -python example.i
gcc -c example_wrap.c example.c  -I$python_inc
# gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace

# tcl 
swig -tcl example.i
gcc -c example_wrap.c example.c 
# gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
gcc -shared -o example.so  example_wrap.o example.o -ltcl -flat_namespace


# swig -c++ -Wall -python H2O.i
# gcc -c H2O_wrap.cxx ../../Library/src/H2O.C ../../Library/src/Fluid.C  -I$python_inc -I../../Library/include
# gcc -shared -o _H2O.so  H2O_wrap.o H2O.o -L$python_lib -lpython3.7m -flat_namespace