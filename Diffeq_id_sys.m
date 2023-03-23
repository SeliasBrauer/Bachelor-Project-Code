function dadt_sum = Diffeq_id_sys(t,x,Xi,nVars,polyorder)
%computes the differenntial equation for a given number of variables and a
%given maximum 
ind = 1;

% poly order 0
g(ind,1) = 1;
ind = ind+1;

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

dadt = (Xi.*g)';
    for i = 1:nVars
        dadt_sum(i,1) = sum(dadt(i,:));
    end
end 
