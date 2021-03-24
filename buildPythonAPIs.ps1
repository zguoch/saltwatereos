# build python api with all python versions in github workflows
# tested for macOS

$basePath_python=$args[0]
$extName_PythonLib=$args[1]
$cmake_config=$args[2]
function buildPythonAPI()
{
    ls
}
# find all python versions
# $basicpath = "C:\Users\zhiku\Downloads\saltwatereos\Library"
foreach($pyversion in Get-ChildItem -Path $basePath_python)
{
    $python_path=$pyversion.FullName
    foreach($pyversion_minor in Get-ChildItem -Path ${python_path}\x64\include)
    {
        $inlucde_path = ${python_path} + "\x64\include\" + $pyversion_minor.BaseName
        $lib_path = ${python_path} + "\x64\libs\"
        echo $inlucde_path
        ls $inlucde_path
        # ls $lib_path
    }
    # FullName
    
}