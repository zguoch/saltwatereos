% 1. build genLUT_2D
% parallel version if you have OpenMP
mex genLUT_2D.cpp -DUSE_OMP=1 -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -I'/usr/local/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl_par -L'/usr/local/lib' -lomp; exit;
% serial version
mex genLUT_2D.cpp -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl; exit;

% 2. build genLUT_3D
% parallel version if you have OpenMP
mex genLUT_3D.cpp -DUSE_OMP=1 -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -I'/usr/local/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl_par -L'/usr/local/lib' -lomp; exit;

% 3. display LUT info
mex infoLUT.cpp -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl; exit;

% 4.1 build lookupLUT_2D
mex lookupLUT_2D.cpp -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl; exit;
% 4.2 build lookupLUT_3D
mex lookupLUT_3D.cpp -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl; exit;
