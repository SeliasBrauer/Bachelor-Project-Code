function Xi = sparsifyDynamics_con(Theta,dXdt,lambda,n,C_in,d)
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

        C = [C; C_vec]; %add constraint vector to constraint matrix
    end 
end
end

k = [];
%{
for i = 1:n/2  
    C_vec = zeros(1, m * n);
    C_vec([ m*((i*2-1) - 1 ) + i*2 , m*(i*2-1) + i*2-1]) = 1;
    k = [k, find(C_vec)];
    %C = [C; C_vec];
    %d = [d; 0 ];
end
%}

% forcing frequency relation
%{
 C_vec = zeros(1, m * n);
 C_vec([k(2),k(3)]) = [0.332/0.1660, -1];
 C = [C;C_vec];
 d = [d;0];

 C_vec = zeros(1, m * n);
 C_vec([k(2),k(6)]) = [3, -1];
 C = [C;C_vec];
 d = [d;0];
%}
 %{
 C_vec = zeros(1, m * n);
 C_vec([k(3)]) = 1;
 C = [C;C_vec];
 d = [d;0];

 C_vec = zeros(1, m * n);
 C_vec([k(4)]) = 1;
 C = [C;C_vec];
 d = [d;0];
 %}
%{
 C_vec = zeros(1, m * n);
 C_vec([k(5)]) = 1;
 C = [C;C_vec];
 d = [d;0];

  C_vec = zeros(1, m * n);
 C_vec([k(6)]) = 1;
 C = [C;C_vec];
 d = [d;0];
%}
 
% stack derivative data
dXdt = reshape(dXdt,[],1);

% compute Sparse regression: sequential least squares
Xi = lsqlin(Theta_copy,dXdt,[],[],C,d,[],[],[],opts);  % initial guess: Least-squares

% lambda is our sparsification knob.
for k=1:10
    
    for i = 1:n
    range = ((i-1)*m+1):( i*m);
    smallinds(range) = abs(Xi(range))< lambda(i); %* mean(abs( nonzeros(Xi(range)) ) ) );   % find small coefficients of individual modes
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
    
    % compute new dynamics with new constraints
    Xi = lsqlin(Theta_copy,dXdt, [],[] ,C,d,[],[],[],opts);
end

%reshape dynamics
Xi = reshape(Xi,[],n);
