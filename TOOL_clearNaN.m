function Matrix = TOOL_clearNaN(Matrix)
    % ȥ��ȫΪNaNֵ����
    rowsWithAllNaN = all(isnan(Matrix), 2);
    Matrix = Matrix(~rowsWithAllNaN, :);
end