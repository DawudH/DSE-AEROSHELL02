classdef aeroProperties
    properties
        alpha;
        cla;
        cda;
        cmya;
    end
    methods
        function obj = aeroProperties()
            A = dlmread('pastille.txt');
            obj.alpha = A(:,1);
            obj.cla = A(:,4);
            obj.cda = A(:,2);
            obj.cmya = A(:,6);
        end
        
        function [cla, cda, cmya] = aeroCoeffs(obj, alpha)
            % [cla, cda, cmya] = aeroCoeffs(obj, alpha) Calculate coeffs
            cla = interp1(obj.alpha, obj.cla, alpha);
            cda = interp1(obj.alpha, obj.cda, alpha);
            cmya = interp1(obj.alpha, obj.cmya, alpha);
        end
    end
end