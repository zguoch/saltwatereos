# build python api with all python versions in github workflows
# tested for macOS

$basePath_python=$args[0]
$HOME_PATH = $args[1]
$BUILD_CONFIGURATION = $args[2]
$BUILD_PLATFORM = $args[3]

echo "============== 传入的参数 ================="
echo $basePath_python
echo $HOME_PATH
echo $BUILD_CONFIGURATION
echo $BUILD_PLATFORM
echo "========================================="
cd $HOME_PATH
cd Library/build

foreach($pyversion in Get-ChildItem -Path $basePath_python)
{
    $python_path = $pyversion.FullName
    $libname = $pyversion.BaseName[0] + $pyversion.BaseName[2]
    $includePath = $python_path + "\include"
    $libPath = $python_path + "\libs\python" + $libname + ".lib"
    $SWIG_EXECUTABLE = $HOME_PATH + "/dependencies_swEOS/windows/swigwin-4.0.2/swig.exe"
    $swigConfig = "-DSWIG_EXECUTABLE=" + $SWIG_EXECUTABLE
    $pyInlcudeConfig = "-DPYTHON_INCLUDE_DIR=" + $includePath
    $pyLibConfig = "-DPYTHON_LIBRARY=" + $libPath
    $buildTypeConfig = "-DCMAKE_BUILD_TYPE=" + $BUILD_CONFIGURATION
    $buildPlatConfig = "-DCMAKE_GENERATOR_PLATFORM=" + $BUILD_PLATFORM
    $buildConfig = "/p:Configuration=" + $BUILD_CONFIGURATION
    echo "=========================== cmake 配置参数 ========================="
    echo $libPath
    echo $SWIG_EXECUTABLE
    echo $swigConfig
    echo $pyInlcudeConfig
    echo $pyLibConfig
    echo $buildTypeConfig
    echo $buildPlatConfig
    echo $buildConfig
    echo "======================================================================"
    # compile
    cmake $swigConfig  $pyInlcudeConfig $pyLibConfig -DBuild_API_tcl=OFF $buildTypeConfig $buildPlatConfig ..
    msbuild /m $buildConfig eosH2ONaCl.vcxproj
    msbuild /m $buildConfig INSTALL.vcxproj

    $pyswEOS_newPath = "../API/python/pyswEOS_python" + $pyversion.BaseName
    mv ../API/python/pyswEOS $pyswEOS_newPath 
    # $statue = "==============================" + $pyversion.BaseName + " done ===================="
    # echo $statue
}

cd ../../