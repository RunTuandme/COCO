function [srp] = light_pressure(x, y, z, r_Sun, CrAM)
% DESCRIPTION:     Calculate the light pressure acceleration.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0

    ae = 6.37813630E+06;
    p0 = 4.5605e-6;         % N/m^2
    AU = 1.495978707e11;    % m
    vector_r = [x; y; z];
    vector_Rs = r_Sun - vector_r;
    
    Rs = norm(vector_Rs);
    r = norm(vector_r);
%     Cr = 1;
%     A = 16;                 % m^2
%     m = 800;                % kg
    delta_r = vector_Rs / Rs;
    cosPhi = dot(vector_r/r, vector_Rs/Rs);
    sinPhi = sqrt(1 - cosPhi^2);
    if ( cosPhi <= 0 ) && ( sinPhi < (ae/r) )
        srp = 0;
    else
        srp = - p0 * (AU / Rs)^2 * CrAM * delta_r;
    end
end