function Xi = sparsifyDynamics_con_mix2(Theta,dXdt,lambda,lambda2,n,C_in)
% faster version for sparsification using constraints. No equality
% constraints are possible using this function. It constrained variables
% will be set equal to zero. 

%takes cell array of constraints C_in
% C_in = {[mode, variable, constant]}

%create copy of candidate library: 
Theta_copy = kron(eye(n), Theta);

m = width(Theta); % number of candindidate functions

smallinds = zeros(1, m * n); %create constraint vector

%unpack constraints C. Create Vector with variable indexes set to small. 
if numel(C_in) >= 1
    for i = 1:numel(C_in)
        con = C_in{i}; %extract value from input cell array 
        var_num = ((con(1)-1)*m + con(2)); %variable number found from mode number and internal variables number
        smallinds(var_num) = true; %give variable factor according to input constraints
    end
end
% Find big indexes
biginds = ~smallinds; 

% stack derivative data
dXdt = reshape(dXdt,[],1);

%create empty coefficient vector
Xi = zeros(m*n,1); 

% Find initial guess using LASSO: 
Xi(biginds) = lasso(Theta_copy(:,biginds),dXdt,"Lambda",lambda);  % initial guess: Least-squares

% lambda is our sparsification knob.
for k=1:10
    
    for i = 1:n
    range = ((i-1)*m+1):( i*m);
    smallinds(range) =  (abs(Xi(range)) < lambda2(i));% * mean(abs( nonzeros(Xi(range)) ) ) );   % find small coefficients of individual modes
    end
    Xi(smallinds == 1) = 0; %threshold small values: 
    %find big indexes
    biginds = ~smallinds; 
    % compute new dynamics with new big index: 
    Xi(biginds) = Theta_copy(:,biginds) \ dXdt ;
end

%reshape dynamics
Xi = reshape(Xi,[],n);
