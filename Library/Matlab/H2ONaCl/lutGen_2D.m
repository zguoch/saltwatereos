% lutGen_2D.m Help file for lutGen_2D MEX file
%
% Generate 2D Adaptive Mesh Refined (AMR) LookUpTable (LUT) in TPX or HPX space.
%   lutGen_2D(xmin, xmax, ymin, ymax, constZ, const_which_var, TorH, outputFile, min_level=4, max_level=6, num_threads=1);
%   * Space is specified by arg TorH: 0 means TPX space; 1 means HPX space.
%   * x/y variable pair is specified by const_which_var, the third variable is constant value specified by constZ, there are three options: 
%     -- (1) 1 means constant T/H, x is X and y is P
%     -- (2) 2 means constant P, x is X and y is T/H
%     -- (3) 3 means constant X, x is T/H and y is P
%   * Output file, the result will be written in .bin file for lookup purpose and .vtu file for visualization using Paraview.
%   * Unit: T[K], P[Pa], X[wt. NacL, 0~1], H[J/kg]
%   * Default args:
%     -- min_level = 4
%     -- max_level = 6
%     -- num_threads = 1, can be set to other int value for parallel computing if the mex file is comparallel compiled
%   * The maximum valid level is 29, but please be careful when set level greater than 12, it will takes much more memory. 
%
%   MEX File function.