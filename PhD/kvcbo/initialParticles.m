function  V = initialParticles(costfunction, KVCBOparam)
            
    N = KVCBOparam.N;
    d = costfunction.d;

    switch KVCBOparam.initData.method
        case 1
            V = randn(d,N);
            V = V./vecnorm(V,2);           
        case 2
            meanDir = KVCBOparam.initData.meanDir;
            k = KVCBOparam.initData.k; 
            V = randVMF(N, meanDir, k)';                   
        case 3     
            V = randn(d,N);
            V = V./vecnorm(V,2);          
            phi = InvHyperSphere(V);
            phi = mod(phi,pi);
            V = HyperSphere(phi);
        otherwise
           % uniform as default
           V = randn(d,N);
           V = V./vecnorm(V,2);     
    end
            
end     

