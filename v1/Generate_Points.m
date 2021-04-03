classdef Generate_Points
    
    methods (Static)
        %%
        function X = NearlyParallel(M, d, ds, Nsp, noise, MOutliers)
            X = zeros(d,Nsp*M+MOutliers);
            U = orth(randn(d,ds));
            epsilon = 0.05;
            
            for i=1:Nsp
                j = randi([1 d],1,1);
                U(j,:) = U(j,:) + epsilon;
                X(:, (i-1)*M+1:i*M) = (U/sqrt(ds)) * randn(ds,M) + sqrt(noise)*randn(d,M)/sqrt(d);
            end
            
            % adding Outliers
            for i=1:MOutliers
                noise = 10;
                X(:,Nsp*M+i) = sqrt(noise)*randn(d,1)/sqrt(d);
            end
        end
        
        %%
        function X = RandomlyChosen(M, d, ds, Nsp, noise, MOutliers)
            X = zeros(d,Nsp*M+MOutliers);
            
            for i=1:Nsp
                U = orth(randn(d,ds));
                X(:, (i-1)*M+1:i*M) = (U/sqrt(ds)) * randn(ds,M) + sqrt(noise)*randn(d,M)/sqrt(d);
            end
            
            % adding Outliers
            for i=1:MOutliers
                noise = 10;
                X(:,Nsp*M+i) = sqrt(noise)*randn(d,1)/sqrt(d);
            end
        end
        
        %%
        function X = LoadDataFromFile(s, field_separator)
            headerLines = 0;
            X = importdata(s, field_separator, headerLines);
        end
        
        %%
        function X = UniformOnSphere(d,M)
           X = randn(d,M);
           X = X./vecnorm(X,2); 
        end
        
        %%
        function X = EquidistantPoints_withOutliers_2D(MList, MOutliers)
            L = size(MList,2);
            Mt = sum(MList);
            X = zeros(2,Mt+MOutliers);

            for i=1:L
                theta_rand = rand(1)*2*pi;
                G_rand = [cos(theta_rand) sin(theta_rand); -sin(theta_rand) cos(theta_rand)];
                A = [linspace(-10, 10, MList(i)); zeros(1,MList(i))];
                X(:, sum(MList(1:i-1))+1:sum(MList(1:i))) = G_rand*A;
            end  
            
            % adding Outliers
            for i=1:MOutliers
                noise = 10;
                X(:,Mt+i) = sqrt(noise)*randn(2,1)/sqrt(2);
            end
        end
        
        %%
        function X = GaussianClusters_2D(MList, mu_List, sigma_List, MOutliers)
            Mt = sum(MList);
            X = zeros(2,Mt);
            L = size(MList,2);
            
            for i=1:L
                xx = mu_List(:,i)*ones(1,MList(i)) + diag(sigma_List(:,i))*randn(2,MList(i));;
                X(:,(i-1)*MList(i)+1:i*MList(i)) = xx;
            end
            
            % adding Outliers
            for i=1:MOutliers
                noise = 10;
                X(:,Mt+i) = sqrt(noise)*randn(2,1)/sqrt(2);
            end
            
        end
        
        %%
        function X = Random_GaussianClusters_2D(MList, MOutliers)
            Mt = sum(MList);
            X = zeros(2,Mt);
            L = size(MList,2);
            mu = 5;         % choose mean uniformly in [-mu, mu]
            sigma = 2;      % choose sigma uniformly in [0,sigma]
            
            for i=1:L
                mu_rand = mu*(1-2*rand(2,1));
                sigma_rand = sigma*rand(2,1);
                xx = mu_rand*ones(1,MList(i)) + diag(sigma_rand)*randn(2,MList(i));
                X(:,(i-1)*MList(i)+1:i*MList(i)) = xx;
            end
            
            % adding Outliers
            for i=1:MOutliers
                noise = 25;
                X(:,Mt+i) = sqrt(noise)*randn(2,1)/sqrt(2);
            end
            
        end

    end
end

