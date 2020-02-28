function makeicon()
{
    mkdir icons.iconset

    sips -z 16 16     swEOS.png --out icons.iconset/icon_16x16.png

    sips -z 32 32     swEOS.png --out icons.iconset/icon_16x16@2x.png

    sips -z 32 32     swEOS.png --out icons.iconset/icon_32x32.png

    sips -z 64 64     swEOS.png --out icons.iconset/icon_32x32@2x.png

    sips -z 64 64     swEOS.png --out icons.iconset/icon_64x64.png

    sips -z 128 128   swEOS.png --out icons.iconset/icon_64x64@2x.png

    sips -z 128 128   swEOS.png --out icons.iconset/icon_128x128.png

    sips -z 256 256   swEOS.png --out icons.iconset/icon_128x128@2x.png

    sips -z 256 256   swEOS.png --out icons.iconset/icon_256x256.png

    sips -z 512 512   swEOS.png --out icons.iconset/icon_256x256@2x.png

    sips -z 512 512   swEOS.png --out icons.iconset/icon_512x512.png

    sips -z 1024 1024   swEOS.png --out icons.iconset/icon_512x512@2x.png

    iconutil -c icns icons.iconset -o logo2.icns
}
function makeicon_doxygon()
{
    mkdir icons.doxygon 
    sips -z 55 55     swEOS_transparent.png --out icons.doxygon/SWEOS55x55.png
    sips -z 1200 1200     swEOS_transparent.png --out icons.doxygon/favicon.png
}
# makeicon

makeicon_doxygon