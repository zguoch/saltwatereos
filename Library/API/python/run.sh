
function distPython()
{
    platTag=$1
    pythonTags=$2
    pkgsPath=$3
    packageName=pyswEOS
    for pythonTag in $pythonTags
    do
        # pythonTag=$(echo $platform | cut -d ":" -f 1)
        # platTag=$(echo $platform | cut -d ":" -f 2)

        majorTag=$(echo $pythonTag | cut -d "." -f 1)
        minorTag=$(echo $pythonTag | cut -d "." -f 2)
        pythonVersion="${pythonTag:0:3}"
        pyTag=py$pythonTag
        echo $pythonTag $pyTag $pythonVersion
        pkgFolder=pyswEOS_python$pythonTag
        cp -rf $pkgsPath/${pkgFolder}* $packageName
        echo $pkgsPath/${pkgFolder}*

        # create __init__.py in module folder !!!important
        cp __init__.py.in ${packageName}/__init__.py

        # include libraries as materials
        # cp MANIFEST.in ${packageName} 

        # change python tag in setup.cfg
        awk -v pythonTag=$pyTag '{gsub(/python-tag = pythonTag/, "python-tag = "pythonTag""); print }' setup.cfg.in  > setup.cfg
        # awk -v pythonVersion=$pythonVersion '{gsub(/python_requires = ==3.6/, "python_requires = =="pythonVersion".*"); print }'
        # build
        python setup.py bdist_wheel --plat-name $platTag

        version=`cat version.txt`
        # upload to PyPi
        whlFile=dist/pyswEOS-$version-$pyTag-none-$platTag.whl
        ls $whlFile
        twine upload $whlFile --verbose

        rm -rf $packageName
        rm -rf ${packageName}.egg-info
    done
}

rm -rf dist build

distPython manylinux2010_x86_64 "2.7 3.5 3.6 3.7 3.8 3.9" /Users/zguo/Downloads/API_Python_Linux
distPython macosx_10_9_x86_64 "2.7 3.5 3.6 3.7 3.8 3.9" /Users/zguo/Downloads/API_Python_MacOSX
# distPython win_amd64 "2.7 3.5 3.6 3.7 3.8 3.9" /Users/zguo/Downloads/API_Python_Windows