
% [needRefine, rho, T] = test_lookupLUT_3D();
% [needRefine, rho, T] = test_lookupLUT_2D();
% test_createLUT_2D();
test_createLUT_3D();

function [needRefine, rho, T] = test_lookupLUT_3D()
    n_sample = 3E3;
    H = (rand(n_sample).*(3.9-0.1) + 0.1).*1E6;
    P = (rand(n_sample).*(2500 - 100) + 100).*1E5;
    X = rand(n_sample).*(1-0.001) + 0.001;
    
    
    [needRefine, rho, T] = lookupLUT_3D('lut_HPX.bin', H, P, X);
    ind = ( needRefine == 1);
   size_sample = size(X);
   num_sample = size_sample(1)*size_sample(2);
   fprintf('%d random points, %d(%.2f%%) points close to phase boundary.\n',num_sample, length(X(ind)), length(X(ind))/num_sample);
end

function [needRefine, rho, T] = test_lookupLUT_2D()
    n_sample = 3E4;
    X = rand(n_sample).*(1-0.001) + 0.001;
    H = (rand(n_sample).*(3.5-0.1) + 0.1).*1E6;
    
    [needRefine, rho, T] = lookupLUT_2D('lut_constP_XH.bin', X, H);
    ind = ( needRefine == 1);
    
%     plot
    plot(X(~ind), H(~ind), 'k.'); hold on;
    plot(X(ind), H(ind), 'ro'); hold off;
end

function test_createLUT_3D()
    num_threads = 8;
    min_level = 2;
    max_level = 8;
    Hmin = 0.1E6; 
    Hmax = 3.9E6;
    Pmin = 100E5;
    Pmax = 2500E5;
    Xmin = 0.001;
    Xmax = 1;
    TorH = 1; % 0 means TPX space; 1 means HPX space
    % test 3D LUT generation
    genLUT_3D(Hmin, Hmax, Pmin, Pmax, Xmin, Xmax, TorH, 'lut_HPX', min_level, max_level, num_threads)
end

function test_createLUT_2D()
    num_threads = 8;
    min_level = 4;
    max_level = 10;
    Hmin = 0.1E6; 
    Hmax = 3.9E6;
    Pmin = 100E5;
    Pmax = 2500E5;
    Xmin = 0.001;
    Xmax = 1;
    TorH = 1; % 0 means TPX space; 1 means HPX space
    % test 2D LUT generation
    genLUT_2D(Xmin, Xmax, Hmin, Hmax, 350E5, 2, TorH, 'lut_constP_XH', min_level, max_level, num_threads)
end



