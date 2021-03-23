# build python api with all python versions in github workflows
# tested for macOS

basePath_python=$1
extName_PythonLib=$2
cmake_config=$3
function buildPythonAPI()
{
    libBuildPath=$1
    includePath=$2
    libPath=$3
    currentPath=${PWD}
    cd $libBuildPath
    cmake cmake -DPYTHON_INCLUDE_DIR=$includePath -DPYTHON_LIBRARY=$libPath $cmake_config ..
    make install
    cd $currentPath
}
# find all python versions
for pyversion in `ls $basePath_python ` 
do
    # echo $pyversion
    python_path=${basePath_python}/${pyversion}
    # find specific version 
    for pyversion_minor in `ls ${python_path}/x64/include`
    do 
        inlucde_path=${python_path}/x64/include/${pyversion_minor}
        lib_path=${python_path}/x64/lib/lib${pyversion_minor}.$extName_PythonLib
        # ls $inlucde_path
        # ls $lib_path
        # build
        buildPythonAPI Library/build $inlucde_path $lib_path
        # mv rename with python version
        mv Library/API/python/pyswEOS Library/API/python/pyswEOS_${pyversion_minor}
        # pring progress 
        echo "==============================" $pyversion_minor " done ===================="
    done 
done