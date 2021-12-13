% 1. build lutGen_2D
% parallel version if you have OpenMP
% mex lutGen_2D.cpp -DUSE_OMP=1 -I'../../../include' -I'/usr/local/include' -L'../../../lib' -leosH2ONaCl_par -L'/usr/local/lib' -lomp; exit;
% serial version
mex lutGen_2D.cpp -I'../../../include' -L'../../../lib' -leosH2ONaCl; exit;

% 2. build lutGen_3D
% parallel version if you have OpenMP
% mex lutGen_3D.cpp -DUSE_OMP=1 -I'../../../include' -I'/usr/local/include' -L'../../../lib' -leosH2ONaCl_par -L'/usr/local/lib' -lomp; exit;
% serial version
mex lutGen_3D.cpp -I'../../../include' -L'../../../lib' -leosH2ONaCl; exit;

% 3. display LUT info
mex lutInfo.cpp -I'../../../include' -L'../../../lib' -leosH2ONaCl; exit;

% 4.1 build lookupLUT_2D
mex lutLookup_2D.cpp -I'../../../include' -L'../../../lib' -leosH2ONaCl; exit;
% 4.2 build lookupLUT_3D
mex lutLookup_3D.cpp -I'../../../include' -L'../../../lib' -leosH2ONaCl; exit;
