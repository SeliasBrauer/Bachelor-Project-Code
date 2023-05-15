function Xi = sparsifyDynamics_con_KKT(Theta,dXdt,lambda,n,C_in,d)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz


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

%find unique constraints from constraint matrix ie. remove repeated constratins: 
[C,ia] = unique(C,'rows'); d = d(ia);  


% stack derivative data
dXdt = reshape(dXdt,[],1);

%implement KKT equation
[KKT_mat_nx ,KKT_mat_ny] = size(Theta_copy'*Theta_copy);

KKT_mat = zeros(KKT_mat_nx + size(C,1), KKT_mat_ny + size(C',2)); 

KKT_mat(1:KKT_mat_nx,1:KKT_mat_ny) = (2*Theta_copy')*Theta_copy;
KKT_mat(KKT_mat_nx+1:end,1:KKT_mat_ny) = C;
KKT_mat(1:KKT_mat_nx,KKT_mat_ny+1:end) = C';

KKT_vec = [(2*Theta_copy')*dXdt; d]; 

% compute Sparse regression: sequential least squares
KKT = KKT_mat \ KKT_vec ; 
Xi = KKT(1:m*n); 
% set small values equal zero: 
Xi(abs(Xi) < 1e-8) = 0; 

% lambda is our sparsification knob.
for k=1:10
    
    for i = 1:n
    range = ((i-1)*m+1):( i*m);
    smallinds(range) =  abs(Xi(range))< lambda(i);% * mean(abs( nonzeros(Xi(range)) ) ) );   % find small coefficients of individual modes
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
    % compute new dynamics with new constraints using KKT: 
    %implement KKT equation
    [KKT_mat_nx ,KKT_mat_ny] = size(Theta_copy'*Theta_copy);

    KKT_mat = zeros(KKT_mat_nx + size(C,1), KKT_mat_ny + size(C',2)); 

    KKT_mat(1:KKT_mat_nx,1:KKT_mat_ny) = (2*Theta_copy')*Theta_copy;
    KKT_mat(KKT_mat_nx+1:end,1:KKT_mat_ny) = C;
    KKT_mat(1:KKT_mat_nx,KKT_mat_ny+1:end) = C';

    KKT_vec = [(2*Theta_copy')*dXdt; d]; 

    % compute New Sparse regression
    KKT = KKT_mat \ KKT_vec ; 
    Xi = KKT(1:m*n); 
    % set small values equal zero: 
    Xi(abs(Xi) < 1e-8) = 0; 
end

%reshape dynamics
Xi = reshape(Xi,[],n);
