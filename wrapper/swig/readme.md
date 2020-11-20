
# SWIG编译python

代码如下

```bash
rm *.o *.so 
rm *_wrap* example.py

python_inc=/Users/zguo/.pyenv/versions/anaconda3-2019.10/include/python3.7m
# 
python_lib=/Users/zguo/.pyenv/versions/anaconda3-2019.10/lib

swig -c -python example.i
gcc -c example_wrap.c example.c  -I$python_inc
# gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
gcc -shared -o _example.so  example_wrap.o example.o -L$python_lib -lpython3.7m -flat_namespace
```
### 注意事项

- 必须指定python动态库的路径和名称，否则会出现编译错误
- 在Mac系统下必须使用`-flat_namespace`参数，否则会出现编译成功但无法在python中用import导入编译得到的python模块
- 最后生成的动态库的名称必须遵守规则！！！`_example.so`前面的下划线不能省略！不能省略！不能省略！[官方的强调如下](http://www.swig.org/Doc1.3/Python.html)：

    When linking the module, **the name of the output file has to match the name of the module prefixed by an underscore**. If the name of your module is "example", then the name of the corresponding object file should be "_example.so" or "_examplemodule.so". The name of the module is specified using the %module directive or the -module command line option.

# Java script

## dependence 

- v8: `brew install v8`

```
# -jsc, -v8 -DV8_VERSION=0x032530
swig -c++ -javascript -v8 example.i
# only once
npm install -g node-gyp

```