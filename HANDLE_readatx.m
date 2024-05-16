function atxdata = HANDLE_readatx(filename)
% DESCRIPTION:     Read the ATX file.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
% NOTES:           GPS:     'G01' - L1                      |            |
%                           'G02' - L2                      |            |
%                           'G05' - L5                      |            |
%                  GLONASS: 'R01' - G1                      |            |
%                           'R02' - G2                      |            |
%                  Galileo: 'E01' - E1                      |            |
%                           'E05' - E5a                     |            |
%                           'E07' - E5b                     |            |
%                           'E08' - E5 (E5a+E5b)            |            |
%                           'E06' - E6                      |            |
%                  Compass: 'C01' - E1                      |            |
%                           'C02' - E2                      |            |
%                           'C07' - E5b                     |            |
%                           'C06' - E6                      |            |
%                  QZSS:    'J01' - L1                      |            |
%                           'J02' - L2                      |            |
%                           'J05' - L5                      |            |
%                           'J06' - LEX                     |            |
%                  SBAS:    'S01' - L1                      |            |
    N_GPS = 32;
    N_BDS = 62;
    N_GAL = 30;
    atxdata = NaN(15, N_GPS + N_BDS + N_GAL);
    fid = fopen(filename, 'r');
    while ~feof(fid)
        line = fgets(fid);
        if contains(line, "START OF ANTENNA")
            while ~contains(line, "END OF ANTENNA")
                line = fgets(fid);
                
                
                if contains(line, "TYPE / SERIAL NO") && line(21) == 'G'
                    prn = str2double(line(22:23));
                end
                if contains(line, "START OF FREQUENCY")
                    if strcmp(line(4:6), 'G01') % L1
                        start_pointer = 1;
                    elseif strcmp(line(4:6), 'G02') % L2
                        start_pointer = 4;
                    elseif strcmp(line(4:6), 'G05') % L5
                        start_pointer = 7;
                    end
                end
                if contains(line, "NORTH / EAST / UP")
                    dcb = str2double({line(4:10) line(14:20) line(23:30)});
                    atxdata(start_pointer: start_pointer+2, prn) = dcb'; 
                    continue
                end
                
                
                if contains(line, "TYPE / SERIAL NO") && line(21) == 'C'
                    prn = str2double(line(22:23)) + N_GPS;
                end
                if contains(line, "START OF FREQUENCY")
                    if strcmp(line(4:6), 'C01') % E1
                        start_pointer = 1;
                    elseif strcmp(line(4:6), 'C02') % E2
                        start_pointer = 4;
                    elseif strcmp(line(4:6), 'C07') % E5b
                        start_pointer = 7;
                    elseif strcmp(line(4:6), 'C06') % E6
                        start_pointer = 10;
                    end
                end
                if contains(line, "NORTH / EAST / UP")
                    dcb = str2double({line(4:10) line(14:20) line(23:30)});
                    atxdata(start_pointer: start_pointer+2, prn) = dcb'; 
                    continue
                end
                
                if contains(line, "TYPE / SERIAL NO") && line(21) == 'E'
                    prn = str2double(line(22:23)) + N_GPS + N_BDS;
                end
                if contains(line, "START OF FREQUENCY")
                    if strcmp(line(4:6), 'E01') % E1
                        start_pointer = 1;
                    elseif strcmp(line(4:6), 'E05') % E5a
                        start_pointer = 4;
                    elseif strcmp(line(4:6), 'E07') % E5b
                        start_pointer = 7;
                    elseif strcmp(line(4:6), 'E08') % E5
                        start_pointer = 10;
                    elseif strcmp(line(4:6), 'E06') % E6
                        start_pointer = 10;
                    end
                end
                if contains(line, "NORTH / EAST / UP")
                    dcb = str2double({line(4:10) line(14:20) line(23:30)});
                    atxdata(start_pointer: start_pointer+2, prn) = dcb'; 
                    continue
                end
            end
        end
    end
end