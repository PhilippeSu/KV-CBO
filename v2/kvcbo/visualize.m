%%
function visualize(V, Va, problem)

    persistent viscounter EE vv Em EM Ne;

    if isempty(viscounter)
        viscounter = 1;
    end

    if viscounter == 1

       figure(1)

       % Deterministic grid for the half-circle
       nn = 6000;
       VV = RandHyperSphere(problem.d,nn,1);

       EE = problem.cost(VV);
       vv = InvHyperSphere(VV);
       Em = min(EE);
       EM = max(EE);
       Ne = EM-Em;
       EE = (EE-Em)/Ne;

       viscounter = viscounter + 1;
    end

    Ep = problem.cost(V);
    Ea = problem.cost(Va);

    Ep = (Ep-Em)/Ne;
    Ea = (Ea-Em)/Ne; 

    v = InvHyperSphere(V);
    va = InvHyperSphere(Va);

    d = problem.d;

    if d == 2

      plot(vv,EE,'.b')
      hold on
      plot(v,Ep,'or', 'MarkerFaceColor', 'r', 'MarkerSize',4)
      hold on
      plot(va, Ea, 'og', 'MarkerFaceColor', 'g', 'MarkerSize',4);
      axis([0 pi -0.2 1.2]);
      hold off

    elseif d == 3

      mesh_res = 40;

      theta = linspace(0,pi,mesh_res);
      phi = linspace(0,pi,mesh_res);
      ve = zeros(d-1,mesh_res*mesh_res);
      for i=1:mesh_res
          ve(1,mesh_res*(i-1)+1:mesh_res*i) = theta(i)*ones(1,mesh_res);
          ve(2,mesh_res*(i-1)+1:mesh_res*i) = phi;
      end

      VE = HyperSphere(ve);
      E = problem.cost(VE);
      E = (E-Em)/Ne;
      Vx = reshape(ve(1,:),[mesh_res,mesh_res]);
      Vy = reshape(ve(2,:),[mesh_res,mesh_res]);
      EE = reshape(E,[mesh_res,mesh_res]);

      subplot(1,2,1)
      surf(Vx,Vy,EE)

      hold on
      scatter3(v(1,:),v(2,:),Ep,'MarkerEdgeColor','k','MarkerFaceColor','r')
      view(2)
      alpha(0.4)
      xlabel('$\theta$')
      ylabel('$\phi$')
      zlabel('$\mathcal{E}_p(\theta,\phi)$')
      axis([0 pi 0 pi -1 1]);
      hold off

      subplot(1,2,2)
      Y = HyperSphere(ve);
      Yx = reshape(Y(1,:),[mesh_res,mesh_res]);
      Yy = reshape(Y(2,:),[mesh_res,mesh_res]);
      Yz = reshape(Y(3,:),[mesh_res,mesh_res]);
      surf(Yx,Yy,Yz,EE)
      alpha(0.4)
      hold on
      scatter3(V(1,:),V(2,:),V(3,:),'MarkerEdgeColor','k','MarkerFaceColor','r')
      xlabel('x')
      ylabel('y')
      zlabel('z')
      axis([-1 1 -1 1 -1 1]);
      drawnow
      hold off

    end

    pause(0.01)

end
