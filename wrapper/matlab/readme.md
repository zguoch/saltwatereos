如果Matlab出现类似下面这样的错误

```
>> mex Fluid.C
Building with 'Xcode with Clang'.
Error using mex
xcodebuild: error: SDK "macosx10.15.4" cannot be located.
xcrun: error: sh -c '/Applications/Xcode.app/Contents/Developer/usr/bin/xcodebuild -sdk macosx10.15.4 -find clang 2> /dev/null'
failed with exit code 16384: (null) (errno=No such file or directory)
xcrun: error: unable to find utility "clang", not a developer tool or in PATH
```

可以通过修改`/Applications/MATLAB_R2018a.app/bin/maci64/mexopts/clang*.xml` 两个文件而排出错误，
```
<SDKVER>
    <cmdReturns name="xcrun -sdk macosx --show-sdk-version | cut -c1-5 "/>
</SDKVER>
```
改为

```
<SDKVER>
    <cmdReturns name="xcrun -sdk macosx --show-sdk-version | cut -c1-5 "/>
</SDKVER>
```