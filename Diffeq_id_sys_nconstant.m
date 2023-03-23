function dadt_sum = Diffeq_id_sys_nconstant(t,x,Xi,nVars,polyorder)
%computes the differenntial equation for a given number of variables and a
%given maximum 
ind = 1;

%skip constant terms
%{
% poly order 0
g(ind,1) = 1;
ind = ind+1;
%}

% poly order 1
for i=1:nVars
    g(ind,1) = x(i);
    ind = ind+1;
end

% poly order 2
if(polyorder>=2)    
    for i=1:nVars
        for j=i:nVars
            g(ind,1) = x(i)*x(j);
            ind = ind+1;
        end
    end
end

% poly order 3
if(polyorder>=3)    
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                g(ind,1) = x(i)*x(j)*x(k);
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
                    g(ind,1) = x(i)*x(j)*x(k)*x(n);
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
                        g(ind,1) = x(i)*x(j)*x(k)*x(n)*x(m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

dadt = (Xi.*g)';
    for i = 1:nVars
        dadt_sum(i,1) = sum(dadt(i,:));
    end
end 
