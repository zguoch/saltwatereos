specie=NaCl
swig -c++ -javascript -node ${specie}.i
rm -rf build
node-gyp configure build
rm ${specie}_wrap.cxx
# test
node test_${specie}.js