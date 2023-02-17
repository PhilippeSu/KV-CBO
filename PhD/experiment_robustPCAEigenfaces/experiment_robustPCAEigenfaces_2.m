%%

clc; clear; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

%% Eigenfaces (with outliers, p = 0.5)

clear; clc; close all;

p = 0.5;

s = '10KUSAdult_Faces_15outliers.txt';

pixelx = 45;                                                % number of pixels in x direcction
pixely = 64;                                                % number of pixels in y direction
d = pixelx*pixely;

if exist(s, 'file') == 0
    disp(' *** File not found ***')
    return
end

X = importdata(s, ',', 0);
X = bsxfun(@minus, X, mean(X,2));                           % subtract the mean

%% SVD

[U,~,~] = svd(X);
postProcess(U(:,1),pixely,pixelx,'Eigenface_SVD_outliers_1.png')
postProcess(-U(:,1),pixely,pixelx,'Eigenface_SVD_outliers_2.png')

%% KV-CBO

costfunction = CostFunction_robustPCAEigenfaces(X,p,0,0,0,d);
params = params_robustPCAEigenfaces(d);

[Va, ~] = KVCBO(costfunction, params);

postProcess(Va,pixely,pixelx,'Eigenface_KVCBO_outliers_1.png')
postProcess(-Va,pixely,pixelx,'Eigenface_KVCBO_outliers_2.png')

%% FMS

if and(exist('fms.m','file'), exist('calc_sdist.m','file'))

    ds=1; % dimension of the subspace to find
    noise=1e-6;
    options=struct;

    options.scaleopt='normal';
    options.p=p;   
    options.initopt='random';
    options.svdopt='normal';
    options.maxiter=100;
    options.epsilon=10^-10;

    L=fms(X',ds,options);

    postProcess(L,pixely,pixelx,'Eigenface_FMS_outliers_1.png')
    postProcess(-L,pixely,pixelx,'Eigenface_FMS_outliers_2.png')

end

%% GD

if isempty(which('spherefactory'))
    error(['Manopt package not found in path.\n' ...
           'Please run importmanopt.']);
       
else
    manifold = spherefactory(d);
    problem.M = manifold;
    problem.cost = @(x) costfunction_ef(X,x,p);
    h = 1e-6;
    problem.egrad = @(x) numericGradient_forward(problem.cost, x, h);

    [x, xcost, info] = steepestdescent(problem); 

    postProcess(x,pixely,pixelx,'Eigenface_GD_outliers_1.png')
    postProcess(-x,pixely,pixelx,'Eigenface_GD_outliers_2.png')

end

%% aux

function postProcess(W,pixely,pixelx,filename)
    P1 = reshape(W, pixely, pixelx);
    minP1 = min(P1(:));
    P1 = double(P1 - minP1) ./ double( max(P1(:)) - minP1 );
    imwrite(im2double(P1), filename);
end

function grad_num = numericGradient_forward(f, x0, h)
  d = size(x0, 1);
  x0 = x0*ones(1,d);
  x0hplus = x0 + h*eye(d);

  grad_num = arrayfun(@(i) (f(x0hplus(:,i)) - f(x0(:,i)))./(h) , 1:d)';
end

function E = costfunction_ef(A,V,p)
     % same function as in robustPCACostFunction.m
    
     delta = 0;
     nn = size(V,2);
     AN = ones(nn,1)*vecnorm(A,2);
     AV = V'*A;
     Z = AN.*AN-AV.*AV;                  % We use the identity||x - <x,v> v ||^2 = ||x||^2 - <x,v>^2
     Z = Z + (delta^2)*ones(size(Z));    % delta = perturbation parameter (for p < 2)
     E = sum((Z)'.^(p/2));
end

