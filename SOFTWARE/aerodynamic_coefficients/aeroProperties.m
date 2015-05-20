classdef aeroProperties
    properties
        alpha;
        cla;
        cda;
        cmya;
    end
    methods
        function obj = aeroProperties(object)
            switch object
                case 'apollo'
                    filestring = 'apollo.txt';
                case 'irve'
                    filestring = 'irve.txt';                    
                case 'isotensoid'
                    filestring = 'isotensoid.txt';
                case 'pastille'
                    filestring = 'pastille.txt';
                case 'torus'
                    filestring = 'torus.txt';    
                case 'ballute'
                    filestring = 'ballute.txt';                    
                otherwise
                    warning(strcat('Case: ', object, ' does not exist'));
                    filestring = 'torus.txt';   
            end
            A = dlmread(filestring);
            obj.alpha = A(:,1);
            obj.cla = A(:,4);
            obj.cda = A(:,2);
            obj.cmya = A(:,6);
        end
        
        function [cla, cda, cmya] = aeroCoeffs(obj, alpha)
%             if sum(alpha > max(obj.alpha))>0 || sum(alpha < min(obj.alpha))>0
%                 error('Outside measured aera!');
%             else
                cla = interp1(obj.alpha, obj.cla, alpha);
                cda = interp1(obj.alpha, obj.cda, alpha);
                cmya = interp1(obj.alpha, obj.cmya, alpha);
%             end
        end
        
        function cla = getCLA(obj, alpha)
            cla = interp1(obj.alpha, obj.cla, alpha);
        end

        function cda = getCDA(obj, alpha)
            cda = interp1(obj.alpha, obj.cda, alpha);
        end
        
        function cmya = getCMYA(obj, alpha)
            cmya = interp1(obj.alpha, obj.cmya, alpha);
        end
        
        function clcd = getCLCD(obj, alpha)
            clcd = obj.getCLA(alpha)./obj.getCDA(alpha);
        end
        
        function cmcl = getCMCL(obj, alpha)
            cmcl = obj.getCMYA(alpha)./obj.getCLA(alpha);
        end
        
        function dCLAdalpha = getLiftGradient(obj, alpha)
            dalpha = 0.1;
            dCLAdalpha = (-obj.getCLA(alpha)+obj.getCLA(alpha+dalpha))/dalpha;
        end
        
        function dCDAdalpha = getDragGradient(obj, alpha)
            dalpha = 0.1;
            dCDAdalpha = (-obj.getCDA(alpha)+obj.getCDA(alpha+dalpha))/dalpha;
        end        
        
        function dCMYAdalpha = getMomentGradient(obj, alpha)
            dalpha = 0.1;
            dCMYAdalpha = (-obj.getCMYA(alpha)+obj.getCMYA(alpha+dalpha))/dalpha;
        end
        
    end
end