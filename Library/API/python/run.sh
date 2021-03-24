
# 拷贝不同系统的Release中的API文件夹到当前目录，分别命名为mac, linux, win
function distPython()
{
    platTag=$1
    pythonTags=$2
    packageName=pyswEOS
    for pythonTag in $pythonTags
    do
        # pythonTag=$(echo $platform | cut -d ":" -f 1)
        # platTag=$(echo $platform | cut -d ":" -f 2)

        majorTag=$(echo $pythonTag | cut -d "." -f 1)
        minorTag=$(echo $pythonTag | cut -d "." -f 2)
        # pyTag=${majorTag}${minorTag}
        pyTag=py$pythonTag
        echo $pythonTag $pyTag
        pkgFolder=pyswEOS_python$pythonTag
        cp -rf $pkgFolder $packageName

        # create __init__.py in module folder !!!important
        touch ${packageName}/__init__.py

        # include libraries as materials
        cp MANIFEST.in ${packageName} 

        # change python tag in setup.cfg
        awk -v pythonTag=$pyTag '{gsub(/python-tag = pythonTag/, "python-tag = "pythonTag""); print }' setup.cfg.in > setup.cfg

        # build
        python setup.py bdist_wheel --plat-name $platTag

        version=`cat version.txt`
        # upload to PyPi
        whlFile=dist/pyswEOS-$version-$pyTag-none-$platTag.whl
        ls $whlFile
        twine upload $whlFile --verbose

        rm -rf $packageName
    done
}

distPython manylinux2010_x86_64 "2.7 3.5m 3.6m 3.7m 3.8 3.9"