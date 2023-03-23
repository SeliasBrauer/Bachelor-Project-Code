function yout = poolData_nconstant(yin,nVars,polyorder)

n = size(yin,1);
ind = 1;

% skip constant terms. 
%{
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;
%}

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

% poly order 2
if(polyorder>=2)    
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

% poly order 3
if(polyorder>=3)    
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

% poly order 4
if(polyorder>=4)    
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for n=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,n);
                    ind = ind+1;
                end
            end
        end
    end
end

% poly order 5
if(polyorder>=5)    
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for n=k:nVars
                    for m=n:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,n).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end
