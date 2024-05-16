function ylist_new = interp_larange_list(tlist, ylist, tlist_new)
% DESCRIPTION:     Lagrange interpolation for a list of data.
% AUTHOR:          ZhangLei
% EMAIL:           the.lei@qq.com
% LAST MODIFIED:   2024-05-15
% VERSION:         1.0
    NT = size(tlist, 1);
    NY = size(ylist, 1);
    NN = length(tlist_new);
    if NY ~= NT
        error('插值序列长度不一致')
    end
    if 11 > NY
        error('插值节点数据不足')
    end
    
    ylist_new = NaN(NN, 1);
    for i = 1: NN
        t = tlist_new(i);
        [index1, index2] = findIndices(tlist, t);
        id_pass = false;
        while (~id_pass)
            try
                start_id = index1 - 4;
                end_id = index2 + 5;
                larange_tlist = tlist(start_id:end_id);
                larange_ylist = ylist(start_id:end_id);
            catch
                if start_id <= 0
                    index1 = index1 + 1;
                    index2 = index2 + 1;
                else
                    index1 = index1 - 1;
                    index2 = index2 - 1;
                end
                continue
            end
            id_pass = true;
        end
        
        ylist_new(i) = lagrange_interp(larange_tlist, larange_ylist, t, 11);
    end
end

function [index1, index2] = findIndices(T, x)
    % 确保T是升序排列的
%     T = sort(T, 'ascend');
    
    % 找到小于等于x的最大值的索引
    index1 = find(T <= x, 1, 'last');
    
    % 如果x小于T中的所有值，则将index1设置为1
    if isempty(index1)
        index1 = 1;
    end
    
    % 找到大于x的最小值的索引
    index2 = find(T > x, 1, 'first');
    
    % 如果x大于T中的所有值，则将index2设置为T的长度
    if isempty(index2)
        index2 = length(T);
    end
    
    % 如果x正好是T中的某个值，则index1和index2会相同，此时需要调整
    if index1 == index2
        if index1 < length(T)
            index2 = index1 + 1;
        else
            index1 = index1 - 1;
        end
    end
end