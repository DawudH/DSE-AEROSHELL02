classdef aeroProperties
    
    properties
        clainterpolant;
        cdainterpolant;
        cmlainterpolant;
        cxainterpolant;
        czainterpolant;
        
        clinterpolant;
        cdinterpolant;
        cmlinterpolant;
        cxinterpolant;
        czinterpolant;        
    end
    methods
        function obj = aeroProperties(name)
            addpath '..\Aero model';
            load('iteration1_1_orbit.mat');
            
            obj.clainterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CRA_aero_array(3,1:end-2), 'linear');
            obj.cdainterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CRA_aero_array(1,1:end-2), 'linear');
            obj.cmlainterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CMA_aero_array(2,1:end-2), 'linear');
            obj.cxainterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CRA_body_array(1,1:end-2), 'linear');
            obj.czainterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CRA_body_array(3,1:end-2), 'linear');

            obj.clinterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CR_aero_array(3,1:end-2), 'linear');
            obj.cdinterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CR_aero_array(1,1:end-2), 'linear');
            obj.cmlinterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CM_aero_array(2,1:end-2), 'linear');
            obj.cxinterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CR_body_array(1,1:end-2), 'linear');
            obj.czinterpolant = griddedInterpolant(mod.alpha_array(1:end-2), mod.CR_body_array(3,1:end-2), 'linear');
        end
        
        function [cl, cd, cml] = aeroCoeffs(obj, alpha)
%             if sum(alpha > max(obj.alpha))>0 || sum(alpha < min(obj.alpha))>0
%                 error('Outside measured aera!');
%             else
                cl = obj.getCL(alpha);
                cd = obj.getCD(alpha);
                cml = obj.getCML(alpha);
%             end
        end
        
        function cla = getCLA(obj, alpha)
            cla = obj.clainterpolant(alpha);
        end

        function cda = getCDA(obj, alpha)
            cda = obj.cdainterpolant(alpha);
        end
        
        function cmya = getCMLA(obj, alpha)
            cmya = obj.cmlainterpolant(alpha);
        end
        
        function cl = getCL(obj, alpha)
            cl = obj.clinterpolant(alpha);
        end

        function cd = getCD(obj, alpha)
            cd = obj.cdinterpolant(alpha);
        end
        
        
        function cd = getCX(obj, alpha)
            cd = obj.cxinterpolant(alpha);
        end
        
        function cd = getCZ(obj, alpha)
            cd = obj.czinterpolant(alpha);
        end
        
        function cd = getCXA(obj, alpha)
            cd = obj.cxainterpolant(alpha);
        end
        
        function cd = getCZA(obj, alpha)
            cd = obj.czainterpolant(alpha);
        end
        
        function cmy = getCML(obj, alpha)
            cmy = obj.cmlinterpolant(alpha);
        end        
        
        function clcd = getCLCD(obj, alpha)
            clcd = obj.getCLA(alpha)./obj.getCDA(alpha);
        end
        
        function cmcl = getCMCL(obj, alpha)
            cmcl = obj.getCMYA(alpha)./obj.getCLA(alpha);
        end
        
        function dCLAdalpha = getLiftGradient(obj, alpha)
            dalpha = 0.01;
            dCLAdalpha = (-obj.getCLA(alpha)+obj.getCLA(alpha+dalpha))/dalpha;
        end
        
        function dCDAdalpha = getDragGradient(obj, alpha)
            dalpha = 0.01;
            dCDAdalpha = (-obj.getCDA(alpha)+obj.getCDA(alpha+dalpha))/dalpha;
        end        
        
        function dCMYAdalpha = getMomentGradient(obj, alpha)
            dalpha = 0.01;
            dCMYAdalpha = (-obj.getCMYA(alpha)+obj.getCMYA(alpha+dalpha))/dalpha;
        end
        
        function cmclcd = getCMCLCD(obj, alpha)
            cmclcd = obj.getCMYA(alpha)./(obj.getCLA(alpha)./obj.getCDA(alpha));
        end
        
    end
end