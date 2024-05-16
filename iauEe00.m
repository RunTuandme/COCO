function ee = iauEe00(date1, date2, epsa, dpsi)
% DESCRIPTION:     Function of SOFA.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
   % Equation of the equinoxes.
   ee = dpsi * cos(epsa) + iauEect00(date1, date2);
end