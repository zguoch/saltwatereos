specie=H2O
# rm example_wrap.cxx
swig -c++ -javascript -node ${specie}.i
rm -rf build
node-gyp configure build
# test
node test_${specie}.js