
[needRefine, rho, T] = test_lookupLUT_3D();
% [needRefine, rho, T] = test_lookupLUT_2D();
% [needRefine, rho, T] = test_lookupLUT_2D_random();
% test_createLUT_2D();
% test_createLUT_3D();

function [needRefine, rho, T] = test_lookupLUT_3D()
    p0 = 200E5;
    x=linspace(0.001, 1, 100);
    h=linspace(0.1, 3.5, 100)*1E6;
    [X,H]=meshgrid(x,h);
    P = H.*0 + p0;
    
    [needRefine, rho, T] = lutLookup_3D('lut_HPX_10.bin', H, P, X);
    % plot
    contourf(X,H,rho);
    % save fig
    saveas(gcf, 'result.pdf');
end

function [needRefine, rho, T] = test_lookupLUT_2D()
    x=linspace(0.001, 1, 100);
    h=linspace(0.1, 3.5, 100)*1E6;
    [X,H]=meshgrid(x,h);
    
    [needRefine, rho, T] = lutLookup_2D('lut_constP_XH_8.bin', X, H);
    % plot
    contourf(X,H,rho);
    % save fig
    saveas(gcf, 'result.pdf');
end

function [needRefine, rho, T] = test_lookupLUT_2D_random()
    n_sample = 500;
    
    X = rand(n_sample).*(1-0.001) + 0.001;
    H = (rand(n_sample).*(3.5-0.1) + 0.1).*1E6;
    
    [needRefine, rho, T] = lutLookup_2D('lut_constP_XH_8.bin', X, H);
    ind = ( needRefine == 1);
%     plot
    plot(X(~ind), H(~ind), 'k.'); hold on;
    plot(X(ind), H(ind), 'ro'); hold off;
end

function test_createLUT_3D()
    num_threads = 8;
    min_level = 2;
    max_level = 6;
    Hmin = 0.1E6; 
    Hmax = 3.9E6;
    Pmin = 100E5;
    Pmax = 2500E5;
    Xmin = 0.001;
    Xmax = 1;
    TorH = 1; % 0 means TPX space; 1 means HPX space
    % test 3D LUT generation
    lutGen_3D(Hmin, Hmax, Pmin, Pmax, Xmin, Xmax, TorH, 'lut_HPX', min_level, max_level, num_threads)
end

function test_createLUT_2D()
    num_threads = 8;
    min_level = 4;
    max_level = 7;
    Hmin = 0.1E6; 
    Hmax = 3.9E6;
    Pmin = 100E5;
    Pmax = 2500E5;
    Xmin = 0.001;
    Xmax = 1;
    TorH = 1; % 0 means TPX space; 1 means HPX space
    % test 2D LUT generation
    lutGen_2D(Xmin, Xmax, Hmin, Hmax, 350E5, 2, TorH, 'lut_constP_XH', min_level, max_level, num_threads)
end



