
rm -rf dist

python setup_mac.py sdist bdist_wheel --plat-name macosx-10.8-x86_64
python setup_win32.py bdist_wheel --plat-name win32
python setup_linux.py sdist bdist_wheel --plat-name manylinux2010_x86_64
# twine upload dist/*.whl --verbose