%%
function Va = computeVa(costfunction, KVCBOparam, V)
    
    % always returns a matrix Va

   global p;

   alpha = KVCBOparam.alpha;
   [d,N] = size(V); 
   
   Es = costfunction.cost(V);
  [Esmin, i] = min(Es);
  
  switch KVCBOparam.logicVa.method
      case 1
          
         Va = sum((V.*exp(-alpha*(Es-Esmin))),2);
         Va = Va./sum(exp(-alpha*(Es-Esmin))); 
         Va = Va*ones(1,N); 
      case 2
          
         Va = V(:,i);
         Va = Va*ones(1,N); 
      case 3
          
         p = KVCBOparam.logicVa.p;
         C = KVCBOparam.logicVa.C; 
         b = KVCBOparam.logicVa.b;

         r = pdist(V', @distfun);
         r = squareform(r);  
         
         Va = zeros(d,N);
         for i = 1:N
            ri = r(i,:);
            omega = 1./(1+ri.^2./C).^(b/2);
            Va(:,i) = sum((V.*exp(-alpha*(Es-Esmin)).*omega),2);
            Va(:,i) = Va(:,i)./sum(exp(-alpha*(Es-Esmin)).*omega); 
         end
         
      case 4
          
          p = KVCBOparam.logicVa.p;
          k = KVCBOparam.logicVa.k;
          r = squareform(pdist(V', @distfun));
          [~,I] = mink(r,k);
          I = I';

          Va = zeros(d,N);
          for i=1:N
              % Idx = knnsearch(V',V(:,i)','K',k,'Distance', @distfun); % 'cosine');
              Eh = Es(I(i,:));
              Ehmin = min(Eh);
              Vh = V(:,I(i,:));
              h = exp(-alpha*(Eh-Ehmin));
              Va(:,i) = sum((Vh.*h),2)/sum(h);
          end
  end
  
end

%% 

function D2 = distfun(ZI,ZJ)
    % ZI is a 1 x n vector
    % ZJ is an m2 x n matrix
    
    global p;
    
    [m2, ~] = size(ZJ);
    D2 = zeros(m2,1);
    for k = 1:m2
       D2(k,1) = norm(ZI - ZJ(k,:) ,p); 
    end
    
end
