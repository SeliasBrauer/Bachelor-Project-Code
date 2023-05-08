function Xi = sparsifyDynamics_con_single(Theta,dXdt,lambda,n,C_in,d)
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

%unpack constraints C
if numel(C_in) >= 1
for i = 1:numel(C_in)
    C_vec = zeros(1, m * n); %create constraint vector
    con = C_in{i}; %extract value from input cell array 
    for j = 1: 3 : length(con) %loop over number of variables being constrained
        var_num = ((con(j)-1)*m + con(j+1)); %variable number found from mode number and internal variables number
        C_vec(var_num) = con(j+2); %give variable factor according to input constraint
    end 
    C = [C; C_vec]; %add constraint vector to constraint matrix
end
end


% stack derivative data
dXdt = reshape(dXdt,[],1);

% compute Sparse regression: sequential least squares
Xi = lsqlin(Theta_copy,dXdt,[],[],C,d,[],[],[],opts);  % initial guess: Least-squares


% lambda is our sparsification knob.
for k=1:50
    smallinds = zeros(length(Xi),1); %vector for storing small indexes. 

    for i = 1:n
    range = ((i-1)*m+1):( i*m); %coefficient range for mode i. 
    Xi(Xi == 0) = NaN;   %temperaly remove zero values: 
    [min_c,min_indx] = min(abs(Xi(range))); % find smallest coefficient in range.
    Xi(isnan(Xi)) = 0; %Change NaN back to zero.  
   

    %if smallest coefficint for mode i is smalled than lambda for mode i.
    %add to small index. 
    smallinds(range(min_indx)) =  min_c< lambda(i); % * mean(abs( nonzeros(Xi(range)) ) ) );   % find small coefficients of individual modes
    end
    
    %smallinds = (abs(Xi)< lambda * mean(abs( nonzeros(Xi) ) ) );   % find small coefficients
    indx = find(smallinds); 

    %add constrainst on small coeffieicnts
    for i = 1:length(indx)
        c_vec = zeros(1,width(Theta_copy)); %create zero vector
        c_vec(indx(i)) = 1; % identify index with small coefficeint
        C = [C; c_vec]; % append to constraint matrix
        d = [d; 0]; % force small coefficient to be zero. 
    end
    
    % compute new dynamics with new constraints
    Xi = lsqlin(Theta_copy,dXdt, [],[] ,C,d,[],[],[],opts);
end

%reshape dynamics
Xi = reshape(Xi,[],n);
