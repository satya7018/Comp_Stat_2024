% Optimal BODs based on Uniform prior

% SEED value
rng( 431)
LB = 0*ones(1,4)';
UB = ones(1,4)';

% Balanced design
xb = [1/4 1/4 1/4 1/4];
x0=xb;
% Initian values to be suppled in the optimization algorithm
x10 = [0.3035    0.2894    0.2635    0.1436];
x20 = [0.2429    0.3288    0.1630    0.2652];
x30 = [0.2298    0.3478    0.1664    0.2560];
x40 = [0.2628    0.3480    0.1564    0.2328];
x50 = [0.2326    0.3213    0.1785    0.2676];

% Parameter space. Use appropriate parameter space. This code uses H3.
lb1 = [-4 0 0 0];
ub1 = [2 3 2 2];

lb2 = [-4 0 0 0 1]; 
ub2 = [2 3 2 2 5];

lb3 = [-4 0 0 0 1.5];
ub3 = [2 3 2 2 3.5];

lb4 = [-4 0 0 0 1];
ub4 = [2 3 2 2  2];

lb5 = [-4 0 0 0 -0.3];
ub5 = [2 3 2 2 0.3];


options = optimoptions('fmincon','TolFun',1e-8,'MaxFunctionEvaluations',10000,'Algorith','interior-point','Display','iter');


% Uncomment to obtain appropriate BOD. This code prodices BOD for Clayton
% Copula

%  [x1,fval1] = fmincon(@(x)BCD_product_Uniform(x,lb1,ub1),x10,[],[],ones(1,4),1,LB,UB,[],options);
 

%  [x2,fval2] = fmincon(@(x)BCD_Frank_Uniform(x,lb2,ub2),x20,[],[],ones(1,4),1,LB,UB,[],options);
 
  % [x3,fval3] = fmincon(@(x)BCD_Gumbel_Uniform(x,lb3,ub3),x30,[],[],ones(1,4),1,LB,UB,[],options)

[x4,fval4] = fmincon(@(x)BCD_Clayton_Uniform(x,lb4,ub4),x40,[],[],ones(1,4),1,LB,UB,[],options)

% [x5,fval5] = fmincon(@(x)BCD_Normal_Uniform(x,lb5,ub5),x50,[],[],ones(1,4),1,LB,UB,[],options);

%  eff1 = (fval1/BCD_product_Uniform(xb,lb1,ub1));

%   eff2 = (fval2/BCD_Frank_Uniform(xb,lb2,ub2));

 % eff3 = (fval3/BCD_Gumbel_Uniform(xb,lb3,ub3))
   eff4 = (fval4/BCD_Clayton_Uniform(xb,lb4,ub4))
% eff5 = (fval5/BCD_Normal_Uniform(xb,lb5,ub5));


% [x1 fval1 eff1; x2 fval2 eff2; x3 fval3 eff3; x4 fval4 eff4; x5 fval5 eff5]




