classdef outputClass
   
    methods (Static)
        
        function [VE, ve] = create_mesh(d) 
            % returns the mesh in Euclidean coordinates (VE)
            % and in spherical coordinates (ve)
            
            mesh_res = 500;
            
            if d == 2
                phi = linspace(0,pi,mesh_res);   
                VE = HyperSphere(phi);
                ve = InvHyperSphere(VE);
                
            elseif d == 3
                theta = linspace(0,pi,mesh_res);
                phi = linspace(0,pi,mesh_res);
                ve = zeros(d-1,mesh_res*mesh_res);
                for i=1:mesh_res
                    ve(1,mesh_res*(i-1)+1:mesh_res*i) = theta(i)*ones(1,mesh_res);
                    ve(2,mesh_res*(i-1)+1:mesh_res*i) = phi;
                end
                
                VE = HyperSphere(ve);
                ve = InvHyperSphere(VE);
            end
            
        end
        
       function visualizeCost(costfunction, VE, ve, V, Va, Vtr)
           
            % rescaling
            EE = costfunction.cost(VE);
            Em = min(EE);
            EM = max(EE);
            Ne = EM-Em;
            EE = (EE-Em)/Ne;
           
            E = costfunction.cost(V);
            Ea = costfunction.cost(Va);
            E = (E-Em)/Ne;
            Ea = (Ea-Em)/Ne; 

            v = InvHyperSphere(V);
            va = InvHyperSphere(Va);

            d = costfunction.d;

            if d == 2
               
              plot(ve,EE,'-b', 'LineWidth',2.0)
              hold on
              plot(v,E,'or', 'MarkerFaceColor', 'r', 'MarkerSize',4)
              hold on
              plot(va, Ea, 'og', 'MarkerFaceColor', 'g', 'MarkerSize',4);
              axis([0 pi -0.2 1.2]);
              
              xticks([0 pi/2 pi])
              xticklabels({'0','\pi/2','\pi'})
              yticks([1/2 1])
              yticklabels({'1/2','1'})
              hold off
              drawnow

              pause(0.1)

            elseif d == 3
                
              mesh_res = sqrt(size(VE,2));
              Vx = reshape(ve(1,:),[mesh_res,mesh_res]);
              Vy = reshape(ve(2,:),[mesh_res,mesh_res]);
              EE = reshape(EE,[mesh_res,mesh_res]);
              
              % First plot
              subplot(1,2,1)
              alpha(1)
              box on
              
              surf(Vx,Vy,EE,'EdgeAlpha','0')
              hold on
              scatter3(v(1,:),v(2,:),E,'MarkerEdgeColor','k','MarkerFaceColor','r')
              view(2)

              xlabel('$\theta$')
              ylabel('$\phi$')
              zlabel('$\mathcal{E}_p(\theta,\phi)$')
              axis([0 pi 0 pi -1 1]);
              hold off

              % Second plot
              subplot(1,2,2)
              alpha(1)
              box on
              
              Y = HyperSphere(ve);
              Yx = reshape(Y(1,:),[mesh_res,mesh_res]);
              Yy = reshape(Y(2,:),[mesh_res,mesh_res]);
              Yz = reshape(Y(3,:),[mesh_res,mesh_res]);

              s = surf(Yx,Yy,Yz,EE, 'EdgeAlpha' , '0');
              direction = [1 -1 0];
              rotate(s,direction,5)     
              hold on
              scatter3(V(1,:),V(2,:),V(3,:),'MarkerEdgeColor','k','MarkerFaceColor','r')
              xlabel('x')
              ylabel('y')
              zlabel('z')
              axis([-1 1 -1 1 0 1.2]);
              box on
              
              %% plot trajectories
              [Ntr,t] = size(Vtr);
              Mtr = 1:d:t;
              Vtr_aux = zeros(d,size(Mtr,2));
              for j = 1:Ntr
                 for l = 0:d-1
                  Vtr_aux(l+1,:) = Vtr(j,Mtr+l); 
                 end                 
                  hold on
                  plot3(Vtr_aux(1,:),Vtr_aux(2,:),Vtr_aux(3,:), 'Color','k', 'LineWidth',1.5)
                  hold on
                  scatter3(Vtr_aux(1,1),Vtr_aux(2,1),Vtr_aux(3,1),'MarkerEdgeColor','k','MarkerFaceColor','y')
              end
              
              %%
              drawnow
              hold off
              pause(0.1)

            end
            
       end
        
       function verbose(info, param, problem, iter)

        nT = param.stopCond.nT;
        N = param.N;

        stats = info(iter);

        cost = stats.cost;
        var = stats.variance;
        con = stats.consensus;
        err = stats.error;
        changeVa = stats.changeVa;
        grnorm = stats.gradnorm;
        stopCond = param.stopCond.name;
        name = problem.name;
        mode = param.mode;
        logicVa = param.logicVa.name;
        init = param.initData.name;

        msg = sprintf([
            ' %s \n', ...
            ' *************************************************\n', ...
            ' Initial data: \t \t %s\n', ...
            ' Noise: \t \t %s\n', ...
            ' Consensus point: \t %s\n', ...
            ' Stopping Condition: \t %s\n', ...
            ' max. Iterations: \t %i\n', ...
            ' Iterations: \t \t %i\n', ...
            ' Number of Particles: \t %i\n', ...
            ' Variance: \t \t %.2e\n', ...
            ' Consensus: \t \t %.2e\n', ...
            ' Change in Va: \t \t %.2e\n', ...
            ' Gradient norm: \t %.2e\n', ...
            ' Cost: \t \t \t %.2e\n', ...
            ' Error: \t \t %.2e\n'], ...
            name, init, mode, logicVa, stopCond, nT, iter, N, var, con, changeVa, grnorm, cost, err);

        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr msg]);
        
    end
    
    end

end