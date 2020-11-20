rm *.o *.so 
rm *_wrap* example.py

python_inc=/Users/zguo/.pyenv/versions/anaconda3-2019.10/include/python3.7m
# 
python_lib=/Users/zguo/.pyenv/versions/anaconda3-2019.10/lib

# swig -python example.i
# gcc -c example_wrap.c example.c  -I$python_inc
# # gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
# gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
