for pyversion_minor in `ls C:/hostedtoolcache/windows/Python`
do 
    echo $pyversion_minor
    cmake --version

    msbuild --version
done