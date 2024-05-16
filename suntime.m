function [lst] = suntime(hour, minut, sec, longtitude)
% DESCRIPTION:     Calculate mean solar time.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0

%     lst=0;
%     if longtitude>180
%         longtitude=longtitude-360;
%     end
%     lst=hour+minut/60+(sec)/3600;
%     lst=lst+(longtitude/360)*24;
    if longtitude > 180
        longtitude = longtitude - 360;
    end
    longtitude_hour_angel = longtitude / 15;
%     longtitude = longtitude - 120;
%     delta_minute = longtitude * 4;
    lst = hour + minut / 60 + sec / 3600 - longtitude_hour_angel;
    if lst < 0
        lst = lst + 24;
    end
end