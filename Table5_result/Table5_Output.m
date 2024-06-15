% This program produce the results for True copula = FC and Parameter space =
% H2, 
rng(464)
mu = [-1 1.5 -1 -1];
sgm = 2;
% Parameter space H2; change the parameter space according to H1, H2, and
% H3. Note that parameter space for alpha is required here.

lb1 = 1; 
ub1 = 3;

lb2 = 1.5;
ub2 = 2.5;

lb3 = 1;
ub3 = 3;

lb4 = -0.3;
ub4 = 0.3;

% Follwoing design are given in Table 8 in Appendix C. Use appropriate
% designs to generate the data in Table 4

% Optimal design based on H3, using uniform prior and Frank copula
x1 = [0.157, 0.186, 0.409, 0.248];

% Optimal design based on H3, using uniform prior and Gumbel copula
x2 = [0.253, 0.075, 0.432, 0.241];

% Optimal design based on H3, using uniform prior and Clayton copula
x3 = [0.233, 0.092, 0.442, 0.233];

% Optimal design based on H3, using uniform prior and Normal copula
x4 = [0.184, 0.180, 0.419, 0.217];


% Optimal design based on H3, using Normal prior and Frank copula
xn1 = [0.128, 0.271, 0.324, 0.277];

% Optimal design based on H3, using  Normal prior and Gumbel copula
xn2 = [0.237, 0.190, 0.375, 0.198];

% Optimal design based on H3, using  Normal prior and Clayton copula
xn3 = [0.213, 0.172, 0.337, 0.278];

% Optimal design based on H3, using  Normal prior and Normal copula
xn4 = [0.189, 0.253, 0.338, 0.220];


fval11 = BCD_Frank_Normal(xn1, mu, lb1, ub1, sgm)/BCD_Frank_Normal(x1, mu, lb1, ub1, sgm);

fval21 = BCD_Gumbel_Normal(xn2, mu, lb2, ub2, sgm)/BCD_Gumbel_Normal(x2, mu, lb2, ub2, sgm)

fval31 = BCD_Clayton_Normal(xn3, mu, lb3, ub3, sgm)/BCD_Clayton_Normal(x3, mu, lb3, ub3, sgm)

fval41 = BCD_Normal_Normal(xn4, mu, lb4, ub4, sgm)/BCD_Normal_Normal(x4, mu, lb4, ub4, sgm);


[fval11 fval21 fval31 fval41]


 







