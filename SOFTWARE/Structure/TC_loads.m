%Input parameters
r_ts = 12;
r_a  = 5 ;
r_t  = 0.5;
theta_tc = 60;
Cd = 1.2;
qmax = 1000;

attachment = 'a'


theta_tc = theta_tc .* pi./180;

Drag_IAD = qmax.*Cd.*pi.*r_ts.^2;


if strcmp(attachment,{'a'})
    p_infl = Drag_IAD./(pi^2*r_t^2*cos(theta_tc));
    sig_ts = [Drag_IAD./(2.*pi.*r_a.*cos(theta_tc)) 0];
    sig_t  = [p_infl.*r_t./2 p_infl.*r_t./2*(2*r_ts-r_t)./(r_ts-r_t)];
elseif strcmp(attachment,{'b'})
    p_infl = 3/2*Drag_IAD./(pi^2*r_t^2*cos(theta_tc));
    sig_ts = [Drag_IAD./(2.*pi.*r_a.*cos(theta_tc)) 0];
    sig_t  = [(p_infl.*pi.^2.*r_t.^2.*cos(theta_tc)+Drag_IAD)/(2.*pi.^2.*r_t.*cos(theta_tc)) p_infl.*r_t./2*(2*r_ts-r_t)./(r_ts-r_t)];
end