% Optimal BODs based on Normal prior

% SEED value

rng(166)

LB = 0*ones(1,4)';
UB = ones(1,4)';
sgm = 3.5;
% Balanced design
xb = [1/4 1/4 1/4 1/4];

x0=xb;
% Initian values to be suppled in the optimization algorithm
x10 = [0.1745    0.2697    0.2735    0.2822];
x20 = [0.1403    0.2982    0.3060    0.2555];
x30 = [0.1745    0.2697    0.2735    0.2822];
x40 = [0.1745    0.2697    0.2735    0.2822];
x50 = [0.1847    0.2778    0.3261    0.2114];

% Parameter space. Use appropriate parameter space. This code uses H3.
mu1 = [-1 1.5 1 1];
mu2 = mu1;
lb2 = 1; 
ub2 = 3;

mu3 = mu1;
lb3 = 1.5;
ub3 = 2.5;

mu4 = mu1;
lb4 =  1;
ub4 =  2;

lb5 =  -0.3;
ub5 =  0.3;
mu5 = mu1;

options = optimoptions('fmincon','TolFun',1e-8,'MaxFunctionEvaluations',10000,'Algorith','interior-point','Display','iter');

% Uncomment to obtain appropriate BOD. This code prodices BOD for Clayton
% Copula
% [x1,fval1] = fmincon(@(x)BCD_product_Normal(x,mu1,sgm),x10,[],[],ones(1,4),1,LB,UB,[],options);
% 
% [x2,fval2] = fmincon(@(x)BCD_Frank_Normal(x,mu2,lb2,ub2,sgm),x20,[],[],ones(1,4),1,LB,UB,[],options);
% 
% [x3,fval3] = fmincon(@(x)BCD_Gumbel_Normal(x,mu3,lb3,ub3,sgm),x30,[],[],ones(1,4),1,LB,UB,[],options);

[x4,fval4] = fmincon(@(x)BCD_Clayton_Normal(x,mu4,lb4,ub4,sgm),x40,[],[],ones(1,4),1,LB,UB,[],options)

% [x5,fval5] = fmincon(@(x)BCD_Normal_Normal(x,mu5,lb5,ub5,sgm),x50,[],[],ones(1,4),1,LB,UB,[],options);

% eff1 = (fval1/BCD_product_Normal(xb,mu1,sgm))
% 
% eff2 = (fval2/BCD_Frank_Normal(xb,mu2,lb2,ub2,sgm))
% 
% eff3 = (fval3/BCD_Gumbel_Normal(xb,mu3, lb3,ub3,sgm))

eff4 = (fval4/BCD_Clayton_Normal(xb, mu4, lb4, ub4, sgm))

% eff5 = (fval5/BCD_Normal_Normal(xb,mu5, lb5,ub5,sgm))






