function Xi = sparsifyDynamics_con(Theta,dXdt,lambda,n,C,d)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

opts = optimset('Display','off');

%create copy of candidate library: 
Theta_copy = kron(eye(n), Theta);

m = width(Theta); % number of candindidate functions
k = [];
for i = 1:n/2  
    C_vec = zeros(1, m * n);
    C_vec([ m*((i*2-1) - 1 ) + i*2 , m*(i*2-1) + i*2-1]) = 1;
    k = [k, find(C_vec)];
    C = [C; C_vec];
    d = [d; 0 ];
end

C_vec = zeros(1, m * n);
C_vec(k) = 1; 
C = [C; C_vec];
d = [d; 0];

% stack derivative data
dXdt = reshape(dXdt,[],1);

% compute Sparse regression: sequential least squares
Xi = lsqlin(Theta_copy,dXdt,[],[],C,d,[],[],[],opts);  % initial guess: Least-squares

% lambda is our sparsification knob.
for k=1:10
    smallinds = (abs(Xi)< lambda * mean(abs( nonzeros(Xi) ) ) );   % find small coefficients
    indx = find(smallinds); 

    %add constrainst on small coeffieicnts
    for i = 1:length(indx)
        c_vec = zeros(1,width(Theta_copy)); %create zero vector
        c_vec(indx(i)) = 1; % identify index with small coefficeint
        C = [C; c_vec]; % append to constraint matrix
        d = [d; 0]; % force coefficient to be zero. 
    end
    
    % compute new dynamics with new constraints
    Xi = lsqlin(Theta_copy,dXdt, [],[] ,C,d,[],[],[],opts);
end

%reshape dynamics
Xi = reshape(Xi,[],n);