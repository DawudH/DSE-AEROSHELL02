function find_orbit(rx,CD,fid,begin,ending,step,refine)
    old_crash = 1; old_atmos = 1; old_orbit = 1; win = 0;
    for i=1:length(rx)
        [R,V,a,conditions] = full_orbit(rx(i),CD);
        disp(['rx = ' num2str(rx(i)) ' [m], CD = ' num2str(CD) ' [-], in atmosphere: ' num2str(conditions.atmos) ', crashed: ' num2str(conditions.crash) ', in orbit: ' num2str(conditions.orbit)])
        fprintf(fid,'%11.1f %4.2f %d %d %d \n',rx(i),CD,conditions.atmos,conditions.crash,conditions.orbit);
        if (conditions.atmos < old_atmos) && (conditions.orbit == false) && (i == 1)
                rx = (2*begin -ending):step:begin;
                find_orbit(rx,CD,fid,begin,ending,step,refine);
        end
        if (conditions.crash < old_crash) && (conditions.orbit == false)
            if i == 1
                rx = (2*begin -ending):step:begin;
                find_orbit(rx,CD,fid,begin,ending,step,refine);
            else
                step = step/refine;
                rx = rx(i-1):step:rx(i);
                find_orbit(rx,CD,fid,begin,ending,step,refine);
            end
        end
        if (conditions.orbit == true)
           win = win + 1;
           if win>10
               break;
           end
        end
        old_crash = conditions.crash; old_atmos = conditions.atmos; old_orbit = conditions.orbit;
    end