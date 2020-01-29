%
% Test Minimum: Example function 
% in arbitrary dimension
%
% computation of the minimum
% by KV-CBO method
%
% 28/01/2020
%
% Compute averages over several runs
%

clear; clc; close all;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 20;  % dimension of the space 

Testf={'Rastrigin function','Ackley function'};

%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%

LOOP = 100;      % Number of averages
STEP = LOOP/10;  % Output intervals
BATCH = 0;       % Choose batch computation

% Minimizers parameters

mu = 0.2 ;   % Fast solver intensity
tf = 10;     % Fast solver intervals

% Numerical parameters

Mp=70;     % points for averaging 
Np=100;    % points in the minimizer
Nmin=10;   % Minimum number of particles

alphav = 50000;
lambda=1;

%%%%%%%%%

answer = questdlg('Please, select function:', ...
	'KV-CBO minimizers', ...
	'Ackley function','Rastrigin function','Exit','Exit');

% Handle response

switch answer
    case 'Ackley function'
        energy_functional = @Energy_Class.Ackley;
        ex = 2; SCHEME=2; sigma=.3; dt=0.05; T=2000;
        disp('**************************************************');
        disp(sprintf('KV-CBO solver: d=%1.0f, N=%1.0f, M=%1.0f, mu=%2.1f, alpha=%1.0f',d,Np,Mp,mu,alphav));
        disp(sprintf('Isotropic scheme: sigma=%2.1f dt=%2.1e It=%1.0f',sigma,dt,T));
        disp('**************************************************');
        
    case 'Rastrigin function'
        energy_functional = @Energy_Class.Rastrigin;
        ex = 1; SCHEME=3; sigma=3; dt=0.0025; T=4000;
        disp('**************************************************');
        disp(sprintf('KV-CBO solver: d=%1.0f, N=%1.0f, M=%1.0f, mu=%2.1f, alpha=%1.0f',d,Np,Mp,mu,alphav));
        disp(sprintf('Componentwise scheme: sigma=%2.1f dt=%2.1e It=%1.0f',sigma,dt,T));
        disp('**************************************************');
        
    case 'Exit'
        return;
end
disp(sprintf('Computing %1.0f averages...',LOOP));

B = [zeros(d-1,1);1];
%B = ones(1,d)'./sqrt(d);
Vmin = B ;
Emin = energy_functional(Vmin,B);

% Loop over averages

succ=zeros(1,LOOP);
Error = 0;
Npavg = 0;
g = zeros(1,T+1);
g2 = zeros(1,T+1);
g3 = zeros(1,T+1);
g4 = zeros(1,T+1);
er = zeros(1,T+1);
nn = zeros(1,T+1);
Np0=Np;


for trun = 1:LOOP

    
    Np=Np0;
    
    if round(trun/STEP)*STEP == trun
        disp(sprintf('Run N.%1.0f',trun));
    end
    
    % Initial data uniform on the half hypersphere
    
    Vp = RandHyperSphere(d,Np,1);
    
    % Random subset
    
    if Mp < Np
        sp = randperm(Np,Mp);
        Vs = Vp(:,sp);
    else
        sp = 1:Np;
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
    
    % Vectorialization of Valpha
    
    Valpha=Valpha*ones(1,Np);
    
    % Errors 

    GG(1)=sum(vecnorm(Vp-Valpha,2))/Np;
    GG2(1)=sum(vecnorm(Vp-Vmin*ones(1,Np),2))/Np;
    GG3(1)=max(abs(Valpha(:,1)-Vmin));
    
    Vbar = mean(Vp');
    GG4(1)=sum(vecnorm(Vp-Vbar'*ones(1,Np),2).^2)/Np;
    
    EER(1) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;
    
    NN(1)=Np;
    
    % KV-CBO minimizer
    
    Vp1 = Vp;
    
    for t=1:T
        
        if BATCH==0
            sp = 1:Np;
        end
        
        Salphi = ones(d,1)*sum(Valpha(:,sp).*Vp(:,sp));
        coll = (Valpha(:,sp)-Salphi.*Vp(:,sp));
        
        Nalphi = ones(d,1)*vecnorm(Vp(:,sp)-Valpha(:,sp),2);
        coll2 = Nalphi.^2.*Vp(:,sp); 
        
        eta = randn(d,Np);
        Setai = ones(d,1)*sum(eta(:,sp).*Vp(:,sp));
        diff = Nalphi.*(eta(:,sp)-Setai.*Vp(:,sp));
        
        coll2 = vecnorm(Vp(:,sp)-Valpha(:,sp),2).^2;
       
        % Main loop
        
        switch SCHEME
            
            case 1 % Euler-Maruyama
                
                Vp1(:,sp) = Vp(:,sp).*(1-dt.*(sigma.^2./2).*coll2.*(d-1)) +...
                    dt*lambda*coll+sqrt(dt).*sigma.*diff;
                
            case 2 % Semi-Implicit scheme
                
                tau = 1./(1+dt.*(sigma.^2./2).*coll2.*(d-1));
                Vp1(:,sp) = tau.*(Vp(:,sp) + dt*lambda*coll+sqrt(dt).*sigma.*diff);
                
            case 3 % Component-wise noise
                
                Nalphi = (Vp(:,sp)-Valpha(:,sp));
                
                coll2 = (Nalphi.^2 + vecnorm(Nalphi,2).^2 - ...
                    2*vecnorm(Nalphi.*Vp(:,sp),2)).*Vp(:,sp);
                
                eta = randn(d,Np);
                Setai = ones(d,1)*sum(Nalphi.*eta(:,sp).*Vp(:,sp));
                diff = Nalphi.*eta(:,sp)-Setai.*Vp(:,sp);
                
                Vp1(:,sp) = Vp(:,sp) + dt*lambda*coll +sqrt(dt).*sigma.*diff-...
                    dt*(sigma.^2./2).*coll2;
        end
        
        % Renormalization on the circle (necessary because of the time discretization error)
        
        NNp = vecnorm(Vp1,2);
        Vp = Vp1./NNp;
        
        % Computing variance
        
        Vbar = mean(Vp')'*ones(1,Np);
        GG4(t+1)=sum(vecnorm(Vp-Vbar,2).^2)/Np;
        
        NN(t+1)=Np;
        
        if (round(t/tf)*tf == t)&&(mu>0)
            
            % Particle reduction
            
            Np1 = max(round(Np*(1+mu*((GG4(t+1)-GG4(t+1-tf))/GG4(t+1-tf)))),Nmin);
            NN(t+1)=Np1;
            
            if Np1 < Np
                Np = Np1;
                Vp = Vp(:,1:Np);
                Vp1 = Vp;
                NN(t+1)=Np;
                
            end
        end
        
        %%%%%%%%
        
        % Random subset
        
        if Mp < Np
            sp = randperm(Np,Mp);
            Vs = Vp(:,sp);
        else
            Vs = Vp;
             sp = 1:Np;
        end
        
        
        % New energy values
        
        Es = energy_functional(Vs,B);
        Ems = min(Es);
        
        % Computation of Valpha
        
        Valpha_old = Valpha;
        
        Valpha = sum((Vs.*exp(-alphav*(Es-Ems)))')';
        Nalpha = sum(exp(-alphav*(Es-Ems)));
        Valpha = Valpha / Nalpha;
        
        % Renormalization on the circle
        
        NNalpha = norm(Valpha,2);
        Valpha = (Valpha/NNalpha);
        
        % Vectorialization of Valpha
        
        Valpha=Valpha*ones(1,Np);
        
        % Errors
        
        GG(t+1)=sum(vecnorm(Vp-Valpha,2))/Np;
        GG2(t+1)=sum(vecnorm(Vp-Vmin*ones(1,Np),2))/Np;
        GG3(t+1)=max(abs(Valpha(:,1)-Vmin));
        EER(t+1) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;
        
    end
    
    GG3(T)=max(abs(Valpha(:,1)-Vmin));
    EER(T) = vecnorm(abs(Valpha(:,1)-Vmin),2).^2;
     
    Npavg = Npavg + mean(NN);
     
    if GG3(T) < 0.25
        succ(trun)=1;
        g = g + GG;
        g2 = g2 + GG2;
        g3 = g3 + GG3;
        g4 = g4 + GG4;
        er = er + EER;
        nn = nn + NN;
        Error = Error + EER(T);
       
    else
        disp(sprintf('No Convergence reached Consensus:%f',GG(T)));
    end
end

Ns = sum(succ);

Error = Error/Ns;
Npavg = Npavg/LOOP ;

disp(sprintf('Success rate: %1.0f/%1.0f',sum(succ),LOOP));
disp(sprintf('Expectation error: %g',Error));
disp(sprintf('Average particles number: %f',Npavg));

GG  = g./Ns;
GG2 = g2./Ns;
GG3 = g3./Ns;
GG4 = g4./Ns;
EER = er./Ns;
NN = nn/Ns;


time = (0:T)*dt;

semilogy(time,GG,'--.','Color',[0.9290 0.6940 0.1250]);
hold on;
semilogy(time,GG2,'--x','Color',[0.4660 0.6740 0.1880]);
semilogy(time,GG3,'--*','Color',[0.8500 0.3250 0.0980]);
semilogy(time,EER,'--o','Color',[0 0.4470 0.7410]);
xlabel('t');
legend('$\frac1{N}\sum_{i=1}^N \|V^i-V^{\alpha,\mathcal{E}}\|_2$','$\frac1{N}\sum_{i=1}^N \|V^i-v^*\|_2$',...
    '$\|V^{\alpha,\mathcal{E}}-v^*\|_\infty$','$\frac1{d}\|V^{\alpha,\mathcal{E}}-v^*\|_2^2$','Location','southwest');
legend('boxoff');
hold off;

box off
