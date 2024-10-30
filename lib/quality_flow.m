function [MAPE,RMSE] = quality_flow(actual,predicted)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
 % 检查输入张量的大小是否相同
    if ~isequal(size(actual), size(predicted))
        error('输入张量的大小必须相同');
    end
    
    % 将三维张量转换为二维矩阵
    actual = reshape(actual, [], size(actual, 3));
    predicted = reshape(predicted, [], size(predicted, 3));
    omega = (actual~=0);
    % 计算百分比误差
    absError = abs(actual.*omega - predicted.*omega);
    percentageError = absError ./ actual;

    % 处理实际值为0的情况
    percentageError(isnan(percentageError)) = 0;
    percentageError(isinf(percentageError)) = 0;

    % 计算MAPE
    MAPE = mean(percentageError(:)) * 100;

    % 计算平方误差
    squaredError = (actual - predicted).^2;
    
    % 计算均方根误差
    RMSE = sqrt(mean(squaredError(:)));
end

