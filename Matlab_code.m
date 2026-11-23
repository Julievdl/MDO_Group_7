% Design vector: z = [z1; z2; .....; zn]
objective = @(z) [];

% Start from the initial values given in assignment   
z0 = [];

% Setting the lower and upper bound
lb = [];
ub = [];

% No linear constraints in this assignment
A = [];  b = [];
Aeq = [];  beq = [];

% Nonlinear constraint function
nonlcon = @aircraftConstraints;

[zopt, fhat_opt, exitflag, output] = fmincon( ...
    objective, z0, A, b, Aeq, beq, lb, ub, nonlcon);

xopt = []; %depends on ratio between z and design variable
Dopt = []; %depends on ratio between fhat_opt and Dopt

function [cineq, ceq] = aircraftConstraints(z)
cineq = [];   % Must be <= 0
ceq = [];      % No equality constraints
end



area_opt = [];
constraint_values = [];

zopt;
fhat_opt;
exitflag;
output.message