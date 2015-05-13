function [ out ] = refine_check(c , rx,in)
%checks if rx needs to be refined

    %checks if you are too close
    if c.crash && c.orbit == false
        to_close = true;
    else
        to_close = false;
    end
    
    %checks if you are too far
    if c.flyby
        to_far = true;
    else
        to_far = false;
    end
    
    %saves the rx at the last time you where too close
    if to_close
        out.rx = rx;
    else
        out.rx = in.rx;
    end
    
    %if you are too far for the first time, refine
    if to_far
        refine = true;
    else
        refine = false;
    end
    
    %%output
    out.refine = refine;
end

