# build python api with all python versions in github workflows
# tested for macOS

$basePath_python=$args[0]
$HOME_PATH = $args[1]
$BUILD_CONFIGURATION = $args[2]
$BUILD_PLATFORM = $args[3]

cd $HOME_PATH
cd Library/build

foreach($pyversion in Get-ChildItem -Path $basePath_python)
{
    $python_path = $pyversion.FullName
    $libname = $pyversion.BaseName[0] + $pyversion.BaseName[2]
    $includePath = $python_path + "\include"
    $libPath = $python_path + "\libs\python" + $libname + ".lib"
    $SWIG_EXECUTABLE = $HOME_PATH + "/dependencies_swEOS/windows/swigwin-4.0.2/swig.exe"
    echo $libPath
    # compile
    cmake -DSWIG_EXECUTABLE=$SWIG_EXECUTABLE -DPYTHON_INCLUDE_DIR=$includePath -DPYTHON_LIBRARY=$libPath -DBuild_API_tcl=OFF -DCMAKE_BUILD_TYPE=$BUILD_CONFIGURATION -DCMAKE_GENERATOR_PLATFORM=$BUILD_PLATFORM ..
    msbuild /m /p:Configuration=$BUILD_CONFIGURATION eosH2ONaCl.vcxproj
    msbuild /m /p:Configuration=$BUILD_CONFIGURATION INSTALL.vcxproj

    $pyswEOS_newPath = "../API/python/pyswEOS_python" + $pyversion.BaseName
    mv ../API/python/pyswEOS $pyswEOS_newPath 

    echo "==============================" $pyversion.BaseName " done ===================="
}

cd ../../