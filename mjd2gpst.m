function [GPSweek, GPSsec] = mjd2gpst(mjd)
% DESCRIPTION:     Convert MJD to GPS week and GPS seconds.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% INPUT:           mjd: Modified Julian Date
    GPSweek = floor((mjd-44244)/7);
    GPSsec = floor((mjd-44244-GPSweek*7)*86400*1000)/1000;
end