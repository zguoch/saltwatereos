
% parallel version if you have OpenMP
mex genLUT_2D.cpp -DUSE_OMP=1 -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -I'/usr/local/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl_par -L'/usr/local/lib' -lomp; exit;

% serial version
mex genLUT_2D.cpp -I'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/include' -L'/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/saltwatereos/Library/lib' -leosH2ONaCl; exit;