options = struct;
options.Display = 'iter'; 
options.maxiter = 200;

% i=1: (20,66), [0.25,2,0.5]
% i=2: (25,88), [1,2,0.5]
% i=3: (18,75), [1,1,1]
% i=4: (20,66), [1,2,0.5]
% i=5; (18,75), [1,1,1]
% i=6; (25,86.5), [1,4,0.25]
fminsearch(@(x) func_3param_345_R47(x,6),ones(3,1),options)