function Xi_big = sparsifyDynamics_con_lasso(Theta,dXdt,lambda,n,C_in,d,alpha)
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
% compute Sparse regression: sequential least squares
Xi = lasso(Theta_copy(:,biginds),dXdt,"Lambda",lambda,"Alpha",alpha);  % initial guess: Least-squares

Xi_big(biginds) = Xi; 

%reshape dynamics
Xi_big = reshape(Xi_big,[],n);
