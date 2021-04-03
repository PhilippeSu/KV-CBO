%
% Test Minimum: Demo file with
% examples function in 3D
%
% computation of the minimum
% by Kuramoto-Vicsek CBO and
% graphic visualization
%
% 25/01/2020
%

clear; clc; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 3;                      % dimension of the space (fixed here)

Testf = {'Rastrigin function','Ackley function'};

%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%

VIDEO = 1;                   % Video frame storage
PLOT = 1;                    % plotting data
RENDER = 1;                  % Rendering of the surface
tras = 0.7;                  % transparency

cp = [0.9290 0.6940 0.1250]; % point color
ca = [0 .5 .8];              % Valpha color
cm = [.8 0 0];               % Minimum color

mesh_resolution = 80;                      % mesh resolution

% Minimizer parameters

SCHEME = 1;                  % Scheme
PERI = 1;                    % Periodization 

% Numerical parameters

lambda = 1;

% Input box for other parameters

prompt = {'Enter number of points N:','Enter batch size M:',...
    'Enter value of \sigma:','Enter value of \Delta t:',...
    'Enter value of \alpha:','Enter number of iterations:','Enter direction of minimum:'};
dlgtitle = 'KV-CBO';
dims = [1 35];
definput = {'30','30','0.7','0.1','500','80','0 0 1'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
if length(answer)==0
    return;
end

N = str2num(answer{1}); 
N3 = str2num(answer{2});
sigma = str2num(answer{3});
dt = str2num(answer{4});
alphav = str2num(answer{5});
T = str2num(answer{6});
B = str2num(answer{7});

B = B'./norm(B,2);

% Input box for examples

answer = questdlg('Please, select function:', ...
	'Functions', ...
	'Ackley function','Rastrigin function','Exit','Exit');
% Handle response
switch answer
    case 'Ackley function'
        energy_functional = @Energy_Class.Ackley;
        ex = 2; zb = 0; zu = 15;
    case 'Rastrigin function'
        energy_functional = @Energy_Class.Rastrigin;
        ex = 1; zb = 0; zu = 40;
    case 'Exit'
        return;
end

% Capturing the screen size

scrsz = get(0,'ScreenSize');

% Storage values

Vps=zeros(T,d,N);
vps=zeros(T,d-1,N);

%%%%%
    
%B = [0;0;1];
%B = [-1/sqrt(3);-1/sqrt(3);1/sqrt(3)];

Vmin = B ;
Emin = energy_functional(Vmin,B);
vmin = InvHyperSphere(Vmin);

if PLOT==1
    
    % Deterministic grid for the half-sphere
    
    theta = linspace(0,pi,mesh_resolution);
    phi = linspace(0,pi,mesh_resolution);
    th=[theta(1)*ones(1,mesh_resolution);phi];
    for i=2:mesh_resolution
        th=[th [theta(i)*ones(1,mesh_resolution);phi]];
    end
    
    % Data for plotting in spherical coordinates
   
    v = th;
    Y = HyperSphere(v);
    
    %x = 6*(v-pi/2)/pi;

    Xx=reshape(v(1,:),[mesh_resolution,mesh_resolution]);
    Xy=reshape(v(2,:),[mesh_resolution,mesh_resolution]);
    
    Yx=reshape(Y(1,:),[mesh_resolution,mesh_resolution]);
    Yy=reshape(Y(2,:),[mesh_resolution,mesh_resolution]);
    Yz=reshape(Y(3,:),[mesh_resolution,mesh_resolution]);
    
    E = energy_functional(Y,B);
    
end

% Initial data uniform on the half hypersphere

Vp = RandHyperSphere(3,N,1);

% Random subset

if N3 < N
    sp = randperm(N,N3);
    Vs = Vp(:,sp);
else
    Vs = Vp;
end

% Initial energy 

Es = energy_functional(Vs,B);

% Actual Energy minimum 

Ems = min(Es);

% Computation of Valpha
    
Valpha = sum((Vs.*exp(-alphav*(Es-Ems)))')';
Nalpha = sum(exp(-alphav*(Es-Ems)));
Valpha = Valpha ./ Nalpha;

% Renormalization on the circle

NNalpha = norm(Valpha,2);
Valpha = (Valpha/NNalpha);

%%%%%

valpha = InvHyperSphere(Valpha);
Ealpha = energy_functional(Valpha,B);

% Vectorialization of Valpha

Valpha=Valpha*ones(1,N);

% Errors 

GG(1)=sum(vecnorm(Vp-Valpha,2))/N;
GG2(1)=sum(vecnorm(Vp-Vmin*ones(1,N),2))/N;
GG3(1)=max(abs(Valpha(:,1)-Vmin));
EER(1) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;

if PLOT==1

    Ep = energy_functional(Vp,B);
    EE=reshape(E,[mesh_resolution,mesh_resolution]);
    
    figure('Position',[scrsz(3)/6 scrsz(4)/6 scrsz(3)/1.5 scrsz(4)/1.5],'Name','KV-CBO method');
    subplot(2,2,1);
    
    vp=InvHyperSphere(Vp);
    surf(Xx,Xy,EE);
    hold on;
    scatter3(vp(1,:),vp(2,:),Ep,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cp);
    scatter3(valpha(1),valpha(2),Ealpha,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',ca);
    scatter3(vmin(1),vmin(2),Emin,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cm);
    xlabel('$\theta$');
    ylabel('$\phi$');
    axis([0 pi 0 pi zb zu])
    alpha(tras);
    colormap hsv;
    
    if RENDER == 1
            
            shading interp
            lightangle(170,30)
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
    end
    view(170,30);
    hold off
    
    subplot(2,2,3);
    contour(Xx,Xy,EE);
    hold on;
    scatter3(vp(1,:),vp(2,:),Ep,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cp);
    scatter3(valpha(1),valpha(2),Ealpha,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',ca);
    scatter3(vmin(1),vmin(2),Emin,72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cm);
    xlabel('$\theta$');
    ylabel('$\phi$');
    axis([0 pi 0 pi])
    view(-180,90);
    hold off
    
    subplot(2,2,[2 4]);
    
    surf(Yx,Yy,Yz,EE);
    alpha(tras);
    colormap hsv;
    
    if RENDER == 1
            
            shading interp
            lightangle(-45,30)
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
    end
    
    axis([-1 1 -1 1 -1 1]);
    axis square
    hold on
    scatter3(Vp(1,:),Vp(2,:),Vp(3,:),72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cp);
    quiver3(0,0,0,Vmin(1,1),Vmin(2,1),Vmin(3,1),0,'-k','LineWidth',3);
    quiver3(0,0,0,Valpha(1,1),Valpha(2,1),Valpha(3,1),0,'-r','LineWidth',3);
    scatter3(Valpha(1,:),Valpha(2,:),Valpha(3,:),72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',ca);
    scatter3(Vmin(1),Vmin(2),Vmin(3),72,'LineWidth',1.5,...
        'MarkerEdgeColor','k','MarkerFaceColor',cm);
    
    xlabel('$v_x$');
    ylabel('$v_y$');
    zlabel('$v_z$');
    title(sprintf('$\\|V^{\\alpha,\\mathcal{E}}-v^*\\|_{\\infty}$=%3.2e',GG3(1)));
    
    hold off;
    sgtitle(sprintf('%s: $v^*=(%2.1f,%2.1f,%2.1f)$, N=%1.0f, n=%1.0f, t=%2.1f',Testf{ex},Vmin(1),Vmin(2),Vmin(3),N,0,0),'FontSize',20);
    drawnow;
end

% Stochastic Kuramoto-Vicksek CBO 

Vps(1,:,:)=Vp;

for t=1:T
    
    if PLOT == 0
        disp(sprintf('It=%1.0f E1=%4.3f E2=%4.3f',t,GG(t),GG2(t)));
    end
    
    Salphi = ones(3,1)*sum(Valpha.*Vp);
    eta = randn(3,N);
    Setai = ones(3,1)*sum(eta.*Vp);
    Nalphi = ones(3,1)*vecnorm(Vp-Valpha,2);  
    
    
    sgn = 1;
    norma = Nalphi;
    
    coll = sgn.*(Valpha-Salphi.*Vp);
    coll2 = vecnorm(Vp-Valpha,2).^2.*Vp./vecnorm(Vp);
    
    diff = norma.*(eta-Setai.*Vp);
    
    % Main loop
        
    switch SCHEME
        
        case 1 % Euler-Maruyama
            Vp1 = Vp + dt*lambda*coll+sqrt(dt)*sigma*diff- ...
                dt*(sigma^2/2)*coll2*(d-1);
            
            % Renormalization on the circle 
            
            NNp = vecnorm(Vp1,2);
            Vp = Vp1./NNp;
            
        case 2 % Semi-Implicit scheme
            coll2 = vecnorm(Vp-Valpha,2).^2;
            tau = 1./(1+dt*(sigma^2/2)*coll2*(d-1));
            Vp1 = tau.*(Vp + dt*lambda*coll+sqrt(dt)*sigma*diff);
            
            % Renormalization on the circle 
            
            NNp = vecnorm(Vp1,2);
            Vp = Vp1./NNp;
        
    end
    
    %%%%%%%%
    
    if PERI
        
        vp=InvHyperSphere(Vp);
        
        % Periodization forces the vector to stay on the halph sphere
        
        tv = vp(2,:)>pi;
        vp(2,:) = (vp(2,:)-pi).*tv+vp(2,:).*(1-tv);
        
        %vp = mod(vp,pi);
        
        Vp = HyperSphere(vp);
    end
    
    % Random subset

    if N3 < N
        sp = randperm(N,N3);
        Vs = Vp(:,sp);
    else
        Vs = Vp;
    end
    
    
    % New energy values

    Es = energy_functional(Vs,B);
    
    Ems = min(Es);
    
    % Computation of Valpha
    
    Valpha = sum((Vs.*exp(-alphav*(Es-Ems)))')';
    Nalpha = sum(exp(-alphav*(Es-Ems)));
    Valpha = Valpha / Nalpha;
    
    
    % Renormalization on the circle
    
    NNalpha = norm(Valpha,2);
    Valpha = (Valpha/NNalpha);
    
    valpha = InvHyperSphere(Valpha);
    Ealpha = energy_functional(Valpha,B);

    % Vectorialization of Valpha
    
    Valpha=Valpha*ones(1,N);
    
    % Errors
    
    GG(t+1)=sum(vecnorm(Vp-Valpha,2))/N;
    GG2(t+1)=sum(vecnorm(Vp-Vmin*ones(1,N),2))/N;
    GG3(t+1)=max(abs(Valpha(:,1)-Vmin));
    EER(t+1) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;
    
    if PLOT==1
        
        Ep = energy_functional(Vp,B);
        EE=reshape(E,[mesh_resolution,mesh_resolution]);
        
        subplot(2,2,1);
        vp=InvHyperSphere(Vp);
        
        surf(Xx,Xy,EE);
        hold on;
        scatter3(vp(1,:),vp(2,:),Ep,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',cp);
        scatter3(valpha(1),valpha(2),Ealpha,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',ca);
        scatter3(vmin(1),vmin(2),Emin,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',cm);
        xlabel('$\theta$');
        ylabel('$\phi$');
        axis([0 pi 0 pi zb zu])
        alpha(tras);
        colormap hsv;
        
        if RENDER == 1
            
            shading interp
            lightangle(170,30)
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
        end
        view(170,30);
    
        hold off
        
        subplot(2,2,3);
        contour(Xx,Xy,EE);
        hold on;
        scatter3(vp(1,:),vp(2,:),Ep,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',cp);
        scatter3(valpha(1),valpha(2),Ealpha,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',ca);
        scatter3(vmin(1),vmin(2),Emin,72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',cm);
        xlabel('$\theta$');
        ylabel('$\phi$');
        axis([0 pi 0 pi])
        view(-180,90);
    
        hold off
        
        subplot(2,2,[2 4]);     
        
        surf(Yx,Yy,Yz,EE);
        alpha(tras);
        
        axis([-1 1 -1 1 -1 1]);
        axis square
        hold on
        scatter3(Vp(1,:),Vp(2,:),Vp(3,:),72,'LineWidth',1.5,...
            'MarkerEdgeColor','k','MarkerFaceColor',cp);
         quiver3(0,0,0,Vmin(1,1),Vmin(2,1),Vmin(3,1),0,'-k','LineWidth',3);
         quiver3(0,0,0,Valpha(1,1),Valpha(2,1),Valpha(3,1),0,'-r','LineWidth',3);
         scatter3(Valpha(1,:),Valpha(2,:),Valpha(3,:),72,'LineWidth',1.5,...
             'MarkerEdgeColor','k','MarkerFaceColor',ca);
         scatter3(Vmin(1),Vmin(2),Vmin(3),72,'LineWidth',1.5,...
             'MarkerEdgeColor','k','MarkerFaceColor',cm);
         xlabel('$v_x$');
        ylabel('$v_y$');
        zlabel('$v_z$');
        
        colormap hsv;
        
        if RENDER == 1
            
            shading interp
            lightangle(-45,30)
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.3;
            h.DiffuseStrength = 0.8;
            h.SpecularStrength = 0.9;
            h.SpecularExponent = 25;
            h.BackFaceLighting = 'unlit';
        end
        title(sprintf('$\\|V^{\\alpha,\\mathcal{E}}-v^*\\|_{\\infty}$=%3.2e',GG3(t+1)));
    
        hold off;
        sgtitle(sprintf('%s: $v^*=(%2.1f,%2.1f,%2.1f)$, N=%1.0f, n=%1.0f, t=%2.1f',...
            Testf{ex},Vmin(1),Vmin(2),Vmin(3),N,t,t*dt),'FontSize',20);
        drawnow;
    end
    % Saving video
    
    if VIDEO
        F(t) = getframe(gcf) ;
    end
    
    Vps(t+1,:,:)=Vp;
    
end

GG3(T)=max(abs(Valpha(:,1)-Vmin));
EER(T) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;
    
if GG3(T) < 0.25
    disp('Convergence reached');
else
    disp('No Convergence reached');
end

time=(0:T)*dt;

figure;
semilogy(time,GG,'--.','Color',[0.9290 0.6940 0.1250]);
hold on;
semilogy(time,GG2,'--x','Color',[0.4660 0.6740 0.1880]);
semilogy(time,GG3,'--*','Color',[0.8500 0.3250 0.0980]);
semilogy(time,EER,'--o','Color',[0 0.4470 0.7410]);
title(sprintf('Convergence behavior for %s',Testf{ex}));
xlabel('t');
legend('$\frac1{N}\sum_{i=1}^N \|V^i-V^{\alpha,\mathcal{E}}\|_2$','$\frac1{N}\sum_{i=1}^N \|V^i-v^*\|_2$',...
    '$\|V^{\alpha,\mathcal{E}}-v^*\|_\infty$','$\frac1{d}\|V^{\alpha,\mathcal{E}}-v^*\|_2^2$','Location','southwest');
legend('boxoff');
hold off;

box off

if VIDEO == 1
    
    answer = questdlg('Save video?', ...
        '', ...
        'Yes','No','No');
    
    switch answer
    case 'Yes'
        writerObj = VideoWriter(sprintf('%s',Testf{ex}),'MPEG-4');
        writerObj.FrameRate = 5;
        %writerObj.Quality = 100;
        
        % set the seconds per image
        % open the video writer
        
        open(writerObj);
        
        % write the frames to the video
        
        for i=1:length(F)
            
            % convert the image to a frame
            
            frame = F(i) ;
            writeVideo(writerObj, frame);
        end
        
        % close the writer object
        
        close(writerObj);
    end
end
