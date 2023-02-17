classdef setUpClass
   
    methods (Static)
        
         function [costfunction, KVCBOparam] = Rastrigin(d)
            [costfunction] = RastriginCostFunction(d);
            KVCBOparam = KVCBOparamForRastrigin(d);
        end
        
        function [costfunction, KVCBOparam] = Ackley(d)
            [costfunction] = AckleyCostFunction(d);
            KVCBOparam = KVCBOparamForAckley(d);
        end
        
        function [costfunction, KVCBOparam] = Alpine(d)
            [costfunction] = AlpineCostFunction(d);
            KVCBOparam = KVCBOparamForAlpine(d);
        end
        
        function [costfunction, KVCBOparam] = Schaffer(d)
            [costfunction] = SchafferCostFunction(d);
            KVCBOparam = KVCBOparamForSchaffer(d);
        end
        
        function [costfunction, KVCBOparam] = Solomon(d)
            [costfunction] = SolomonCostFunction(d);
            KVCBOparam = KVCBOparamForSolomon(d);
        end
        
        function [costfunction, KVCBOparam] = Levi(d)
            [costfunction] = LeviCostFunction(d);
            KVCBOparam = KVCBOparamForLevi(d);
        end
        
        function [costfunction, KVCBOparam] = XSYrandom(d)
            [costfunction] = XSYrandomCostFunction(d);
            KVCBOparam = KVCBOparamForXSYrandom(d);
        end
        
        function [costfunction, KVCBOparam] = Griewank(d)
            [costfunction] = GriewankCostFunction(d);
            KVCBOparam = KVCBOparamForGriewank(d);
        end
       
    end
end

