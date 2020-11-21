# rm example_wrap.cxx
swig -c++ -javascript -node H2O.i
rm -rf build
node-gyp configure build
# node-gyp configure build 