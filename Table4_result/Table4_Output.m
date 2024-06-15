% This program produce the results for True copula = FC and Parameter space =
% H3, 
rng(67)
% Parameter space H3; change the parameter space according to H1, H2, and
% H3

lb1 = [-4 0 0 0 1]; 
ub1 = [2 3 2 2  5];

lb2 = [-4 0 0 0 1.5];
ub2 = [2 3 2 2  3.5];

lb3 = [-4 0 0 0 1];
ub3 = [2 3 2 2  3];

lb4 = [-4 0 0 0 -0.3];
ub4 = [2 3 2 2  0.3];

% Follwoing design are given in Table 8 in Appendix C. Use appropriate
% designs to generate the data in Table 4

% Optimal design based on H3, using uniform prior and Frank copula
x1 = [0.235, 0.308, 0.171, 0.286];

% Optimal design based on H3, using uniform prior and Gumbel copula
x2 = [0.218, 0.344, 0.168, 0.270];

% Optimal design based on H3, using uniform prior and Clayton copula
x3 = [0.280, 0.332, 0.141, 0.247];

% Optimal design based on H3, using uniform prior and Normal copula
x4 = [0.233, 0.321, 0.178, 0.268];



fval11 = BCD_Frank_Uniform(x1,lb1,ub1)/BCD_Frank_Uniform(x1,lb1,ub1);

fval21 = BCD_Gumbel_Uniform(x2,lb2,ub2)/BCD_Gumbel_Uniform(x1,lb2,ub2);

fval31 = BCD_Clayton_Uniform(x3,lb3,ub3)/BCD_Clayton_Uniform(x1,lb3,ub3);

fval41 = BCD_Normal_Uniform(x4,lb4,ub4)/BCD_Normal_Uniform(x1,lb4,ub4);



A = [fval11 fval21 fval31 fval41]


 







