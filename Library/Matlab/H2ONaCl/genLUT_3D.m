% genLUT_3D.m Help file for genLUT_3D MEX file
%
% Generate 3D Adaptive Mesh Refined (AMR) LookUpTable (LUT) in TPX or HPX space.
%   genLUT_3D(xmin, xmax, ymin, ymax, constZ, const_which_var, TorH, outputFile, min_level=4, max_level=6, num_threads=1);
%   * Space is specified by arg TorH: 0 means TPX space; 1 means HPX space.
%   * Please not that x, y, z order MUST BE T/H, P, X
%   * Output file, the result will be written in .bin file for lookup purpose and .vtu file for visualization using Paraview.
%   * Unit: T[K], P[Pa], X[wt. NacL, 0~1], H[J/kg]
%   * Default args:
%     -- min_level = 4
%     -- max_level = 6
%     -- num_threads = 1, can be set to other int value for parallel computing if the mex file is comparallel compiled
%   * The maximum valid level is 29, but please be careful when set level greater than 12, it will takes much more memory. 
%
%   MEX File function.