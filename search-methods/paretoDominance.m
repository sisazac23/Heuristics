function [ND] = paretoDominance(F)
% Simone (2020). Pareto filtering (https://www.mathworks.com/matlabcentral/fileexchange/50477-pareto-filtering), MATLAB Central File
[i, dim] = size(F);
idxs = [1 : i]';
while i >= 1
    old_size = size(F,1);
    indices = sum( bsxfun( @ge, F(i,:), F ), 2 ) == dim;
    indices(i) = false;
    F(indices,:) = [];
    idxs(indices) = [];
    i = i - 1 - (old_size - size(F,1)) + sum(indices(i:end));
end 
ND = F;
end