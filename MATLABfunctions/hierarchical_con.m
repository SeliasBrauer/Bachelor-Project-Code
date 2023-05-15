function [C,d] = hierarchical_con(nVars, polyorder)
%Create hierarchial constraints for system with nVars variables; 
% with polyorder.

% create matrix of individual coefficient locations 
for i = 1:nVars
    for j = 1:nVars + 1 
        if j <= 1 
            L(j,i) = i + (nVars+2-j)*(j-1) ;
        elseif j <= 1 + i 
            L(j,i) = L(j-1,i) + (nVars+2-j);
        else
            L(j,i) = L(j-1,i) +1; 
        end
    end
end

% build constraint cell array: 
C = {}; 
m = height(L);
for i = 1: nVars
    if i == 1 || i == 2
        for j = 3:nVars
            for k = 1:m
                C{end+1} = [i,L(k,j),1];
            end
        end
    elseif mod(i,2) == 1 %is odd
        for j = i+1:nVars
            for k = 1:m
                C{end+1} = [i,L(k,j),1];
            end
        end
    else
        range = [i-1,i+1:nVars];
        for j = range
            for k = 1:m
                C{end+1} = [i,L(k,j),1];
            end
        end
    end
end
d = zeros(size(C,2),1);


end