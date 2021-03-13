for platform in mac:macosx-10.9-x86_64 linux:manylinux2010_x86_64 win32:win32
do
    platName=$(echo $platform | cut -d ":" -f 1)
    platTag=$(echo $platform | cut -d ":" -f 2)
    # create __init__.py in module folder !!!important
	touch ${platName}/pyswEOS/__init__.py
    # include libraries as materials
    cp MANIFEST.in ${platName} 
    # cp versions and setup.py
    cp version.txt ${platName} 
    cp setup.py ${platName} 
    # generate wheel binary package
    cd ${platName} 
    version=`cat version.txt`
    python setup.py sdist bdist_wheel --plat-name $platTag
    twine upload dist/*${version}*.whl --verbose
    cd ..
done
