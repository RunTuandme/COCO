function [al, alm] = gravityfield(phi, lam, x, y, z, order, dc, ds, dc0, JD1, JD2)
% DESCRIPTION:     Calculate the acceleration of gravity field.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           JD1, the large part of the Julian Date.
%                  JD2, the small part of the Julian Date.
%                  x, y, z, satellite position.
%                  order, order of the gravitational field model.
%                  phi, lam, the geocentric latitude and longitude.
%                  dc, ds, dc0, coefficients of tide model.
    global clm0
    global clm
    global slm
    global GM
    global ae
    global gravity_type

    if gravity_type == 3
        global tc0 tc ts acosc0 acosc acoss asinc0 asinc asins
        mjd = JD1 - 2400000.5 + JD2;
        mjd0 = 55197;
        T = 365.25;
        dt = (mjd - mjd0) / T;
        clm0cp = clm0 + tc0 * dt + acosc0 * cos(2*pi*dt) + asinc0 * sin(2*pi*dt) + dc0';
        clmcp = clm + tc * dt + acosc * cos(2*pi*dt) + asinc * sin(2*pi*dt) + dc;
        slmcp = slm + ts * dt + acoss * cos(2*pi*dt) + asins * sin(2*pi*dt) + ds;
    else
        clm0cp = clm0 + dc0';
        clmcp = clm + dc;
        slmcp = slm + ds;
    end
    
    % 计算al
    sphi = sin(phi);
    R = [x; y; z];
    r = norm(R);
    ae_r = ae / r;
    [pl, dpl] = average_Pl_u_mex(sphi, order);
        
    %%
    A = clm0cp(2:order) .* ( ae_r.^(2:order) ) / r^3;
    B = (3:order+1)'.* pl(3:order+1)' + sphi * dpl(2:order)';
    M = A * B;
    N = A * r * dpl(2:order)';
    al = -(M * R - N * [0;0;1]) * GM;
       
    %%
%     sum = 0;
%     for l = 2:70
%         p1 = clm0(l) * (ae/r)^l * (1/r)^3;
%         p2 = (l+1) * pl(l+1);
%         p3 = sphi * dpl(l);
%         p4 = r * dpl(l) * [0;0;1];
%         sum = sum + p1 * ( (p2+p3)*R - p4 );
%     end
%     al = -sum * GM;
   %% 
    % 计算alm
    G = [-sin(lam); cos(lam); 0];
    k = [0; 0; 1];
    [p, p_] = average_Plm_u_mex(sphi, order);
    [cosm, sinm] = cos_sin_m(lam, order);
    %%
    ncosm = repmat(cosm, order-1, 1);
    nsinm = repmat(sinm, order-1, 1);
    MM = clmcp(2:order,1:order) .* ncosm + slmcp(2:order,1:order) .* nsinm;
    NN = clmcp(2:order,1:order) .* nsinm - slmcp(2:order,1:order) .* ncosm;
    JA = ae_r.^meshgrid((2:order),(1:order))' / r^3 .* MM;
    J = sum(sum( JA .*(meshgrid((3:order+1),(1:order))' .* p(2:order,1:order) + sphi * p_(2:order,1:order)) ));
    K = sum(sum( JA * (-r) .* p_(2:order,1:order) ));
    D = sum(sum( meshgrid((1:order),(2:order)) / norm([x;y]) .* ae_r.^meshgrid((2:order),(1:order))' / r .* p(2:order,1:order) .* NN ));
    alm = -(J * R + K * k + D * G) * GM;
    
    %%
%     s = 0;
%     for l = 2:70
%         for m = 1:l
%             p1 = (ae/r)^(l)*(1/r)^3;
%             p2 = (l+1) * p(l,m);
%             p3 = sphi * p_(l,m);
%             p4 = r * p_(l,m);
%             p5 = clm(l, m) * cosm(m) + slm(l, m) * sinm(m);
%             p6 = m * (ae/r)^(l) / norm([x;y]) / r * p(l,m);
%             p7 = clm(l, m) * sinm(m) - slm(l, m) * cosm(m);
% 
%             part1 = (p2+p3) * R - p4 * k;
%             part2 = p5;
%             part3 = p6 * p7 * G;
% 
%             s = s + (p1 * part1 * part2 + part3);
%         end
%     end
%     alm = -s * GM;
%     
end

% % 计算pl,dpl
% function [pl, dpl] = average_Pl_u(u, order)
%     pl = zeros(1, order+1);
%     dpl= zeros(1, order);
%     pl(1) = 1;
%     pl(2) = 1.732050807568877 * u;
%     for l = 2:order
%         part1 = sqrt( (2*l+1)/(2*l-1) );
%         part2 = (2-1/l) * u * pl(l);
%         part3 = sqrt( (2*l-1)/(2*l-3) );
%         part4 = (1-1/l) * pl(l-1);
%         pl(l+1) = part1 * (part2 - part3 * part4);
%     end
%     dpl(1) = 1.732050807568877;
%     for l = 2:order
%         pp1 = l / (1-u^2);
%         pp2 = sqrt( (2*l+1)/(2*l-1) ) * pl(l);
%         pp3 = u * pl(l+1);
%         dpl(l) = pp1 * (pp2 - pp3);
%     end
% end
% 
% % 计算plm，dplm
% function [p, p_] = average_Plm_u(u, order)
%     p = zeros(order+1,order+1);
%     p_ = zeros(order,order);
%     
%     p(1,1) = 1.732050807568877 * sqrt(1-u^2);
%     for l = 2:order+1
%         p(l,l) = sqrt( (2*l+1)/(2*l) ) * sqrt(1-u^2) * p(l-1,l-1);
%     end
%     for l=2:order+1
%         p(l,l-1) = sqrt((2*l+1)) * u * p(l-1,l-1);
%     end
%     for l = 3:order+1
%         for m = 1:l-2
%             p(l,m) = sqrt((2*l+1)*(2*l-1)/(l+m)/(l-m)) * u * p(l-1,m) -...
%                 sqrt((2*l+1)*(l-1+m)*(l-1-m)/(2*l-3)/(l+m)/(l-m)) * p(l-2,m);
%         end
%     end
%     
%     for l = 1:order
%         for m = 1:l
%             p_(l,m) = 1 / sqrt(1-u^2)...
%                 * (sqrt( (l+m+1)*(l-m) )...
%                 * p(l,m+1) - m * u / sqrt(1-u^2) * p(l,m));
%         end
%     end
% end
