function Matrix = TOOL_clearNaN(Matrix)
    % 去除全为NaN值的行
    rowsWithAllNaN = all(isnan(Matrix), 2);
    Matrix = Matrix(~rowsWithAllNaN, :);
end