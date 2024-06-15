theta1 = [1 0.3 0.5 -0.5];
theta2 = [1 0.3 0.5 -0.5 0.5];
theta3 = [1 0.3 0.5 -0.5 1.5];
theta4 = [1 0.3 0.5 -0.5 0.5];
theta5 = [0.2 0.3 0.1 -0.5 0.2];

LB1 =[-inf -inf -inf -inf];
UB1 = [inf inf inf inf];

LB2 =[-inf -inf -inf -inf -1];
UB2 = [inf inf inf inf 3];

LB3 =[-inf -inf -inf -inf -1];
UB3 = [inf inf inf inf 3];

LB4 =[-inf -inf -inf -inf -1];
UB4 = [inf inf inf inf 3];

LB5 =[-inf -inf -inf -inf -1];
UB5 = [inf inf inf inf 1];

% Define the path of directory appropriately

T = readmatrix('/Users/satya/Desktop/COST_Reproducible_Codes/Figure_1/Likelihood_data.xlsx');
options = optimoptions('fmincon','TolFun',1e-18,'MaxFunctionEvaluations',10000,'Algorith','interior-point');
x1 = [];
x2 = [];
x3 = [];
x4 = [];
x5 = [];
T(:, 1) = [];
for i = 1:10000
[x1(i,:),fval1(i)] = fmincon(@(theta)likelihood_Product(theta,T(i,:)),theta1,[],[],[],[],LB1,UB1,[],options);
[x2(i,:),fval2(i)] = fmincon(@(theta)likelihood_Frank(theta,T(i,:)),theta2,[],[],[],[],LB2,UB2,[],options);
[x3(i,:),fval3(i)] = fmincon(@(theta)likelihood_Gumbel(theta,T(i,:)),theta3,[],[],[],[],LB3,UB3,[],options);
[x4(i,:),fval4(i)] = fmincon(@(theta)likelihood_Clyton(theta,T(i,:)),theta4,[],[],[],[],LB4,UB4,[],options);
[x5(i,:),fval5(i)] = fmincon(@(theta)likelihood_Normal(theta,T(i,:)),theta5,[],[],[],[],LB5,UB5,[],options);
end
% Estimates 
% Product Copula
mean(x1)
% Other Copulas
[mean(x2);mean(x3);mean(x4);mean(x5)]



