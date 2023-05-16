function Xi = sparsifyDynamics_con_mix(Theta,dXdt,lambda,lambda2,n,C_in,d,alpha)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

opts = optimset('Display','off');

%takes cell array of constraints C_in
% C_in = {[mode, variable, constant, mode, variable, constant , ...]  [...]  []  }


%create copy of candidate library: 
Theta_copy = kron(eye(n), Theta);

m = width(Theta); % number of candindidate functions

C = [];  %create  empty matrxi for storing constraint vectors
var_num = []; % create empty vector for storing constraned variable numbers
%unpack constraints C
if numel(C_in) >= 1
for i = 1:numel(C_in)
    C_vec = zeros(1, m * n); %create constraint vector
    con = C_in{i}; %extract value from input cell array 
    for j = 1: 3 : length(con) %loop over number of variables being constrained
        var_num(i) = ((con(j)-1)*m + con(j+1)); %variable number found from mode number and internal variables number
        C_vec(var_num(i)) = con(j+2); %give variable factor according to input constraint
    end 
    C = [C; C_vec]; %add constraint vector to constraint matrix
end
end

%find unique constraints from constraint matrix ie. remove repeated constratins: 
[C,ia] = unique(C,'rows'); d = d(ia);  

% stack derivative data
dXdt = reshape(dXdt,[],1);

Xi_big = zeros(m*n, 1);

biginds = ones(m*n, 1);

biginds(var_num) = 0; 
biginds = biginds > 0;

% compute Sparse regression: use lasso for initial guess
Xi = lasso(Theta_copy(:,biginds),dXdt,"Lambda",lambda,"Alpha",alpha);  % initial guess: Least-squares
Xi_big(biginds) = Xi; 

Xi = Xi_big; 

% lambda is our sparsification knob.
for k=1:10
    
    for i = 1:n
    range = ((i-1)*m+1):( i*m);
    smallinds(range) =  abs(Xi(range)) <= lambda2(i); % * mean(abs( nonzeros(Xi(range)) ) ) );   % find small coefficients of individual modes
    end
    
    %smallinds = (abs(Xi)< lambda * mean(abs( nonzeros(Xi) ) ) );   % find small coefficients
    indx = find(smallinds); 

    %add constrainst on small coeffieicnts
    for i = 1:length(indx)
        c_vec = zeros(1,width(Theta_copy)); %create zero vector
        c_vec(indx(i)) = 1; % identify index with small coefficeint
        C = [C; c_vec]; % append to constraint matrix
        d = [d; 0]; % force coefficient to be zero. 
    end

    %find unique constraints from constraint matrix ie. remove repeated constratins: 
    [C,ia] = unique(C,'rows'); d = d(ia);  
    
    % compute new dynamics with new constraints
    Xi = lsqlin(Theta_copy,dXdt, [],[] ,C,d,[],[],[],opts);
end


%reshape dynamics
Xi = reshape(Xi,[],n);
