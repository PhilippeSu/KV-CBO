%% KV-CBO method for robust subspace detection energy 
%  and visualization (only for dimension d = 3)

clear; clc; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

% Parameters for Visualization
getVideo = 1;
mesh_resolution = 30; 
tras = 0.4;             % transparency

% Parameters for KV-CBO method
d = 3;                  % dimension (fixed)
alpha_KV = 1e5;        % alpha parameter
sigma = 0.8;
lambda = 1;
dt = 0.25;
N = 10;                     % number of particles
N2 = N;                     % number of particles to be evolved (randomly chosen)
N3 = N;                     % number of particles to compute Valpha (randomly chosen)

% Parameters for Energy
energy_functional = @Energy_Class.RobustSubspaceDetection;  
p = 2;
delta = 0;

% Parameters for Point Cloud
M = 100; 
ds = 1; 
Nsp = 5; 
noise = 0.01;
MOutliers = 0;

% Generate Point Cloud
% X = Generate_Points.NearlyParallel(M, d, ds, Nsp, noise, MOutliers);
X = Generate_Points.RandomlyChosen(M, d, ds, Nsp, noise, MOutliers);
X = bsxfun(@minus, X, mean(X,2));

% Compute Minimum by SVD
if p == 2
    [Vg,D] = eig(X*X');
    [~,i] = max(diag(D));
    Vmin = Vg(:,i);
end

%% KV-CBO method 

% Initial particles
V = RandHyperSphere(d,N,1);
v = InvHyperSphere(V);
v = mod(v,pi);
Valpha = ComputeValpha(energy_functional, X, V, N, alpha_KV, p, delta);
Ep = energy_functional(X,V,p,delta);

% Energy rescaling
Em = min(Ep);
EMM = max(Ep);
Ne = EMM-Em;
Ep = (Ep-Em)/Ne;

% Visualization
theta = linspace(0,pi,mesh_resolution);
phi = linspace(0,pi,mesh_resolution);
ve = zeros(d-1,mesh_resolution*mesh_resolution);
for i=1:mesh_resolution
    ve(1,mesh_resolution*(i-1)+1:mesh_resolution*i) = theta(i)*ones(1,mesh_resolution);
    ve(2,mesh_resolution*(i-1)+1:mesh_resolution*i) = phi;
end
 
VE = HyperSphere(ve);
E = energy_functional(X,VE,p,delta);
E = (E-Em)/Ne;
Vx = reshape(ve(1,:),[mesh_resolution,mesh_resolution]);
Vy = reshape(ve(2,:),[mesh_resolution,mesh_resolution]);
EE = reshape(E,[mesh_resolution,mesh_resolution]);

figure('Renderer', 'painters', 'Position', [10 10 900 300])

subplot(1,2,1)
surf(Vx,Vy,EE)
alpha(tras)
Ax = axis;
% axis([0 pi 0 pi Ax(5) Ax(6)]);
hold on
scatter3(v(1,:),v(2,:),Ep,'MarkerEdgeColor','k','MarkerFaceColor','r')
view(2)
xlabel('$\theta$')
ylabel('$\phi$')
zlabel('$\mathcal{E}_p(\theta,\phi)$')

subplot(1,2,2)
Y = HyperSphere(ve);
Yx = reshape(Y(1,:),[mesh_resolution,mesh_resolution]);
Yy = reshape(Y(2,:),[mesh_resolution,mesh_resolution]);
Yz = reshape(Y(3,:),[mesh_resolution,mesh_resolution]);
surf(Yx,Yy,Yz,EE)
alpha(tras)
hold on
scatter3(V(1,:),V(2,:),V(3,:),'MarkerEdgeColor','k','MarkerFaceColor','r')
xlabel('x')
ylabel('y')
zlabel('z')

disp('Press any key to start')
pause
     
t = 1;
T = 50;
while (t < T)
    
    [Valpha, V] = KuramotoIteration(energy_functional, d, Valpha, V, X, ...
            alpha_KV, sigma, lambda, dt, N, N2, N3, ...
            p, delta); 
        
    v = InvHyperSphere(V);
    v = mod(v,pi);
    Ep = energy_functional(X,V,p,delta);
    Ep = (Ep-Em)/Ne; 
        
    va = InvHyperSphere(Valpha);
    Ea = energy_functional(X,Valpha,p,delta);  
    Ea = (Ea-Em)/Ne;
        
    subplot(1,2,1)
    surf(Vx,Vy,EE)
    title(sprintf('Time Steps: %i', t))
    alpha(tras)
    Ax = axis;
    % axis([0 pi 0 pi Ax(5) Ax(6)]);
    hold on
    scatter3(v(1,:),v(2,:),Ep,'MarkerEdgeColor','k','MarkerFaceColor','r')
    hold on
    scatter3(va(1,:),va(2,:),Ea,'MarkerEdgeColor','k','MarkerFaceColor','g')
    view(2)
    xlabel('$\theta$')
    ylabel('$\phi$')
    zlabel('$\mathcal{E}_p(\theta,\phi)$')
    drawnow
    hold off
    
    subplot(1,2,2)
    surf(Yx,Yy,Yz,EE)
    alpha(tras)
    hold on
    scatter3(V(1,:),V(2,:),V(3,:),'MarkerEdgeColor','k','MarkerFaceColor','r')
    hold on
    scatter3(Valpha(1,:),Valpha(2,:),Valpha(3,:),'MarkerEdgeColor','k','MarkerFaceColor','g')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
    axis([-1 1 -1 1 -1 1]);
    axis square
    drawnow
    hold off
    
    % pause(0.5)
    drawnow
    if getVideo == 1
        F(t) = getframe(gcf);
    end
    t = t + 1;
end

if p == 2
   error = min([norm(Valpha-Vmin,2), norm(Valpha+Vmin,2)]) 
end

if getVideo == 1
     
    writerObj = VideoWriter(sprintf('%s','Video_'),'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    
    for i=1:length(F)
        frame = F(i) ;
        writeVideo(writerObj, frame);
    end
    
    close(writerObj);
end

