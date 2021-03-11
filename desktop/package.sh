macdeployqt build/swEOS.app -dmg


# 1. sign app
# codesign --deep -fs "Developer ID Application: Zhikui GUO (AFCJBNSQCL)" --options=runtime --timestamp swEOS.app

# # 2. check sign
# codesign -vvv --deep --strict swEOS.app

# 3. generate zip file
# 第三步在finder中直接用mac自带的zip压缩，再到第四步上传，命令行制作的zip似乎用不了
# 

# 4. upload
# xcrun altool --notarize-app -t osx --primary-bundle-id="gzk.swEOS" -u "zhikuiguo@live.cn" -p "app-specific-password" -f swEOS.zip --verbose

# # 5. check status
# xcrun altool --notarization-info 2723c659-11fd-4599-87ce-8585e73f13ce -u "zhikuiguo@live.cn" -p "app-specific-password"

# 6. Staple the notarization ticket
# xcrun stapler staple swEOS.app

# 7. 
# spctl -a -v swEOS.app

# 签名-公证完成后，可以用命令生成dmg文件
# hdiutil create -srcfolder swEOS.app swEOS.dmg