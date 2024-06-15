function y = Clayton_Inf_Matrix(theta,x)
%theta = [1 2 0.5 0.4 2];
%x = [0.3 0.2 0.4 0.1];
eta1_AA = theta(1) + theta(2) + theta(3);
eta2_AA = theta(1) - theta(2) + theta(3) + theta(4);

eta1_AB = theta(1) + theta(2) + theta(3);
eta2_AB = theta(1) - theta(2) - theta(3) + theta(4);

eta1_BA = theta(1) + theta(2) - theta(3);
eta2_BA = theta(1) - theta(2) + theta(3) - theta(4);

eta1_BB = theta(1) + theta(2) - theta(3);
eta2_BB = theta(1) - theta(2) - theta(3) - theta(4);

pi1_AA = exp(eta1_AA)/(1+exp(eta1_AA));
pi2_AA = exp(eta2_AA)/(1+exp(eta2_AA));

pi1_AB = exp(eta1_AB)/(1+exp(eta1_AB));
pi2_AB = exp(eta2_AB)/(1+exp(eta2_AB));

pi1_BA = exp(eta1_BA)/(1+exp(eta1_BA));
pi2_BA = exp(eta2_BA)/(1+exp(eta2_BA));

pi1_BB = exp(eta1_BB)/(1+exp(eta1_BB));
pi2_BB = exp(eta2_BB)/(1+exp(eta2_BB));

% Partial derivatives of pi1's wrt \nu
Dpi1_t1_AA =  pi1_AA*(1-pi1_AA);
Dpi1_t1_AB =  pi1_AB*(1-pi1_AB);
Dpi1_t1_BA =  pi1_BA*(1-pi1_BA);
Dpi1_t1_BB =  pi1_BB*(1-pi1_BB);

% Partial derivatives of pi2's wrt \nu
Dpi2_t1_AA =  pi2_AA*(1-pi2_AA);
Dpi2_t1_AB =  pi2_AB*(1-pi2_AB);
Dpi2_t1_BA =  pi2_BA*(1-pi2_BA);
Dpi2_t1_BB =  pi2_BB*(1-pi2_BB);


% Partial derivatives of pi1's wrt \beta^{*}
Dpi1_t2_AA =  pi1_AA*(1-pi1_AA);
Dpi1_t2_AB =  pi1_AB*(1-pi1_AB);
Dpi1_t2_BA =  pi1_BA*(1-pi1_BA);
Dpi1_t2_BB =  pi1_BB*(1-pi1_BB);

% Partial derivatives of pi2's wrt \beta^{*}
Dpi2_t2_AA =  -pi2_AA*(1-pi2_AA);
Dpi2_t2_AB =  -pi2_AB*(1-pi2_AB);
Dpi2_t2_BA =  -pi2_BA*(1-pi2_BA);
Dpi2_t2_BB =  -pi2_BB*(1-pi2_BB);


% Partial derivatives of pi1's wrt \tau^{*}
Dpi1_t3_AA =  pi1_AA*(1-pi1_AA);
Dpi1_t3_AB =  pi1_AB*(1-pi1_AB);
Dpi1_t3_BA =  -pi1_BA*(1-pi1_BA);
Dpi1_t3_BB =  -pi1_BB*(1-pi1_BB);

% Partial derivatives of pi2's wrt \tatu^{*}
Dpi2_t3_AA =  pi2_AA*(1-pi2_AA);
Dpi2_t3_AB =  -pi2_AB*(1-pi2_AB);
Dpi2_t3_BA =  pi2_BA*(1-pi2_BA);
Dpi2_t3_BB =  -pi2_BB*(1-pi2_BB);

% Partial derivatives of pi1's wrt \gamma^{*}
Dpi1_t4_AA =  0;
Dpi1_t4_AB =  0;
Dpi1_t4_BA =  0;
Dpi1_t4_BB =  0;

% Partial derivatives of pi2's wrt \gamma^{*}
Dpi2_t4_AA =  pi2_AA*(1-pi2_AA);
Dpi2_t4_AB =  pi2_AB*(1-pi2_AB);
Dpi2_t4_BA =  -pi2_BA*(1-pi2_BA);
Dpi2_t4_BB =  -pi2_BB*(1-pi2_BB);


R1_AA = ((pi1_AA)^(-theta(5)) + (pi2_AA)^(-theta(5))-1)^((-1/theta(5)) - 1);
R1_AB = ((pi1_AB)^(-theta(5)) + (pi2_AB)^(-theta(5))-1)^((-1/theta(5)) - 1);
R1_BA = ((pi1_BA)^(-theta(5)) + (pi2_BA)^(-theta(5))-1)^((-1/theta(5)) - 1);
R1_BB = ((pi1_BB)^(-theta(5)) + (pi2_BB)^(-theta(5))-1)^((-1/theta(5)) - 1);

R2_AA_t1 = ((pi1_AA)^(-theta(5)-1))*Dpi1_t1_AA + ((pi2_AA)^(-theta(5)-1))*Dpi2_t1_AA;
R2_AB_t1 = ((pi1_AB)^(-theta(5)-1))*Dpi1_t1_AB + ((pi2_AB)^(-theta(5)-1))*Dpi2_t1_AB;
R2_BA_t1 = ((pi1_BA)^(-theta(5)-1))*Dpi1_t1_BA + ((pi2_BA)^(-theta(5)-1))*Dpi2_t1_BA;
R2_BB_t1 = ((pi1_BB)^(-theta(5)-1))*Dpi1_t1_BB + ((pi2_BB)^(-theta(5)-1))*Dpi2_t1_BB;
% Derivative of p_11 wrt \nu
D_p11_AA_t1 = R1_AA*R2_AA_t1;
D_p11_AB_t1 = R1_AB*R2_AB_t1;
D_p11_BA_t1 = R1_BA*R2_BA_t1;
D_p11_BB_t1 = R1_BB*R2_BB_t1;



R2_AA_t2 = ((pi1_AA)^(-theta(5)-1))*Dpi1_t2_AA + ((pi2_AA)^(-theta(5)-1))*Dpi2_t2_AA;
R2_AB_t2 = ((pi1_AB)^(-theta(5)-1))*Dpi1_t2_AB + ((pi2_AB)^(-theta(5)-1))*Dpi2_t2_AB;
R2_BA_t2 = ((pi1_BA)^(-theta(5)-1))*Dpi1_t2_BA + ((pi2_BA)^(-theta(5)-1))*Dpi2_t2_BA;
R2_BB_t2 = ((pi1_BB)^(-theta(5)-1))*Dpi1_t2_BB + ((pi2_BB)^(-theta(5)-1))*Dpi2_t2_BB;
% Derivative of p_11 wrt \beta^{*}
D_p11_AA_t2 = R1_AA*R2_AA_t2;
D_p11_AB_t2 = R1_AB*R2_AB_t2;
D_p11_BA_t2 = R1_BA*R2_BA_t2;
D_p11_BB_t2 = R1_BB*R2_BB_t2;

R2_AA_t3 = ((pi1_AA)^(-theta(5)-1))*Dpi1_t3_AA + ((pi2_AA)^(-theta(5)-1))*Dpi2_t3_AA;
R2_AB_t3 = ((pi1_AB)^(-theta(5)-1))*Dpi1_t3_AB + ((pi2_AB)^(-theta(5)-1))*Dpi2_t3_AB;
R2_BA_t3 = ((pi1_BA)^(-theta(5)-1))*Dpi1_t3_BA + ((pi2_BA)^(-theta(5)-1))*Dpi2_t3_BA;
R2_BB_t3 = ((pi1_BB)^(-theta(5)-1))*Dpi1_t3_BB + ((pi2_BB)^(-theta(5)-1))*Dpi2_t3_BB;
% Derivative of p_11 wrt \tau^{*}
D_p11_AA_t3 = R1_AA*R2_AA_t3;
D_p11_AB_t3 = R1_AB*R2_AB_t3;
D_p11_BA_t3 = R1_BA*R2_BA_t3;
D_p11_BB_t3 = R1_BB*R2_BB_t3;

R2_AA_t4 = ((pi1_AA)^(-theta(5)-1))*Dpi1_t4_AA + ((pi2_AA)^(-theta(5)-1))*Dpi2_t4_AA;
R2_AB_t4 = ((pi1_AB)^(-theta(5)-1))*Dpi1_t4_AB + ((pi2_AB)^(-theta(5)-1))*Dpi2_t4_AB;
R2_BA_t4 = ((pi1_BA)^(-theta(5)-1))*Dpi1_t4_BA + ((pi2_BA)^(-theta(5)-1))*Dpi2_t4_BA;
R2_BB_t4 = ((pi1_BB)^(-theta(5)-1))*Dpi1_t4_BB + ((pi2_BB)^(-theta(5)-1))*Dpi2_t4_BB;
% Derivative of p_11 wrt \gamma^{*}
D_p11_AA_t4 = R1_AA*R2_AA_t4;
D_p11_AB_t4 = R1_AB*R2_AB_t4;
D_p11_BA_t4 = R1_BA*R2_BA_t4;
D_p11_BB_t4 = R1_BB*R2_BB_t4;

% Derivative of p_10 wrt \theta
D_p10_AA_t1 = Dpi1_t1_AA - D_p11_AA_t1;
D_p10_AB_t1 = Dpi1_t1_AB - D_p11_AB_t1;
D_p10_BA_t1 = Dpi1_t1_BA - D_p11_BA_t1;
D_p10_BB_t1 = Dpi1_t1_BB - D_p11_BB_t1;

D_p10_AA_t2 = Dpi1_t2_AA - D_p11_AA_t2;
D_p10_AB_t2 = Dpi1_t2_AB - D_p11_AB_t2;
D_p10_BA_t2 = Dpi1_t2_BA - D_p11_BA_t2;
D_p10_BB_t2 = Dpi1_t2_BB - D_p11_BB_t2;

D_p10_AA_t3 = Dpi1_t3_AA - D_p11_AA_t3;
D_p10_AB_t3 = Dpi1_t3_AB - D_p11_AB_t3;
D_p10_BA_t3 = Dpi1_t3_BA - D_p11_BA_t3;
D_p10_BB_t3 = Dpi1_t3_BB - D_p11_BB_t3;

D_p10_AA_t4 = Dpi1_t4_AA - D_p11_AA_t4;
D_p10_AB_t4 = Dpi1_t4_AB - D_p11_AB_t4;
D_p10_BA_t4 = Dpi1_t4_BA - D_p11_BA_t4;
D_p10_BB_t4 = Dpi1_t4_BB - D_p11_BB_t4;

% Derivative of p_01 wrt \theta
D_p01_AA_t1 = Dpi2_t1_AA  - D_p11_AA_t1;
D_p01_AB_t1 = Dpi2_t1_AB - D_p11_AB_t1;
D_p01_BA_t1 = Dpi2_t1_BA - D_p11_BA_t1;
D_p01_BB_t1 = Dpi2_t1_BB - D_p11_BB_t1;

D_p01_AA_t2 = Dpi2_t2_AA - D_p11_AA_t2;
D_p01_AB_t2 = Dpi2_t2_AB - D_p11_AB_t2;
D_p01_BA_t2 = Dpi2_t2_BA - D_p11_BA_t2;
D_p01_BB_t2 = Dpi2_t2_BB - D_p11_BB_t2;

D_p01_AA_t3 = Dpi2_t3_AA - D_p11_AA_t3;
D_p01_AB_t3 = Dpi2_t3_AB - D_p11_AB_t3;
D_p01_BA_t3 = Dpi2_t3_BA - D_p11_BA_t3;
D_p01_BB_t3 = Dpi2_t3_BB - D_p11_BB_t3;

D_p01_AA_t4 = Dpi2_t4_AA - D_p11_AA_t4;
D_p01_AB_t4 = Dpi2_t4_AB - D_p11_AB_t4;
D_p01_BA_t4 = Dpi2_t4_BA - D_p11_BA_t4;
D_p01_BB_t4 = Dpi2_t4_BB - D_p11_BB_t4;

% computations for derivatives wrt alpha

C_AA = ((pi1_AA)^(-theta(5)) + (pi2_AA)^(-theta(5)) - 1)^(-1/theta(5));
C_AB = ((pi1_AB)^(-theta(5)) + (pi2_AB)^(-theta(5)) - 1)^(-1/theta(5));
C_BA = ((pi1_BA)^(-theta(5)) + (pi2_BA)^(-theta(5)) - 1)^(-1/theta(5));
C_BB = ((pi1_BB)^(-theta(5)) + (pi2_BB)^(-theta(5)) - 1)^(-1/theta(5));

a_AA = ((pi1_AA)^(-theta(5)))*log(pi1_AA) + ((pi2_AA)^(-theta(5)))*log(pi2_AA);
a_AB = ((pi1_AB)^(-theta(5)))*log(pi1_AB) + ((pi2_AB)^(-theta(5)))*log(pi2_AB);
a_BA = ((pi1_BA)^(-theta(5)))*log(pi1_BA) + ((pi2_BA)^(-theta(5)))*log(pi2_BA);
a_BB = ((pi1_BB)^(-theta(5)))*log(pi1_BB) + ((pi2_BB)^(-theta(5)))*log(pi2_BB);

D_C_AA = (C_AA/theta(5))*(((C_AA)^(theta(5)))*a_AA - log(C_AA));
D_C_AB = (C_AB/theta(5))*(((C_AB)^(theta(5)))*a_AB - log(C_AB));
D_C_BA = (C_BA/theta(5))*(((C_BA)^(theta(5)))*a_BA - log(C_BA));
D_C_BB = (C_BB/theta(5))*(((C_BB)^(theta(5)))*a_BB - log(C_BB));

D_p11_AA_t5 = D_C_AA;
D_p11_AB_t5 = D_C_AB;
D_p11_BA_t5 = D_C_BA;
D_p11_BB_t5 = D_C_BB;

D_p10_AA_t5 = -D_p11_AA_t5;
D_p10_AB_t5 = -D_p11_AB_t5;
D_p10_BA_t5 = -D_p11_BA_t5;
D_p10_BB_t5 = -D_p11_BB_t5;

D_p01_AA_t5 = -D_p11_AA_t5;
D_p01_AB_t5 = -D_p11_AB_t5;
D_p01_BA_t5 = -D_p11_BA_t5;
D_p01_BB_t5 = -D_p11_BB_t5;


p11_AA = C_AA;
p10_AA = pi1_AA - p11_AA;
p01_AA = pi2_AA - p11_AA;

p11_AB = C_AB;
p10_AB = pi1_AB - p11_AB;
p01_AB = pi2_AB - p11_AB;

p11_BA = C_BA;
p10_BA = pi1_BA - p11_BA;
p01_BA = pi2_BA - p11_BA;

p11_BB = C_BB;
p10_BB = pi1_BB - p11_BB;
p01_BB = pi2_BB - p11_BB;


DG_P_AA = diag([1/p11_AA 1/p10_AA 1/p01_AA]);
DG_P_AB = diag([1/p11_AB 1/p10_AB 1/p01_AB]);
DG_P_BA = diag([1/p11_BA 1/p10_BA 1/p01_BA]);
DG_P_BB = diag([1/p11_BB 1/p10_BB 1/p01_BB]);

J_AA = (1/(1-p11_AA-p10_AA-p01_AA))*ones(3);
J_AB = (1/(1-p11_AB-p10_AB-p01_AB))*ones(3);
J_BA = (1/(1-p11_BA-p10_BA-p01_BA))*ones(3);
J_BB = (1/(1-p11_BB-p10_BB-p01_BB))*ones(3);

I_AA = DG_P_AA + J_AA;
I_AB = DG_P_AB + J_AB;
I_BA = DG_P_BA + J_BA;
I_BB = DG_P_BB + J_BB;

DP_AA = [D_p11_AA_t1 D_p10_AA_t1 D_p01_AA_t1; D_p11_AA_t2 D_p10_AA_t2 D_p01_AA_t2;...
    D_p11_AA_t3 D_p10_AA_t3 D_p01_AA_t3; D_p11_AA_t4 D_p10_AA_t4 D_p01_AA_t4;...
    D_p11_AA_t5 D_p10_AA_t5 D_p01_AA_t5];

DP_AB = [D_p11_AB_t1 D_p10_AB_t1 D_p01_AB_t1; D_p11_AB_t2 D_p10_AB_t2 D_p01_AB_t2;...
    D_p11_AB_t3 D_p10_AB_t3 D_p01_AB_t3; D_p11_AB_t4 D_p10_AB_t4 D_p01_AB_t4;...
    D_p11_AB_t5 D_p10_AB_t5 D_p01_AB_t5];

DP_BA = [D_p11_BA_t1 D_p10_BA_t1 D_p01_BA_t1; D_p11_BA_t2 D_p10_BA_t2 D_p01_BA_t2;...
    D_p11_BA_t3 D_p10_BA_t3 D_p01_BA_t3; D_p11_BA_t4 D_p10_BA_t4 D_p01_BA_t4;...
    D_p11_BA_t5 D_p10_BA_t5 D_p01_BA_t5];

DP_BB = [D_p11_BB_t1 D_p10_BB_t1 D_p01_BB_t1; D_p11_BB_t2 D_p10_BB_t2 D_p01_BB_t2;...
    D_p11_BB_t3 D_p10_BB_t3 D_p01_BB_t3; D_p11_BB_t4 D_p10_BB_t4 D_p01_BB_t4;...
    D_p11_BB_t5 D_p10_BB_t5 D_p01_BB_t5];

M_AA = DP_AA*I_AA*DP_AA';
M_AB = DP_AB*I_AB*DP_AB';
M_BA = DP_BA*I_BA*DP_BA';
M_BB = DP_BB*I_BB*DP_BB';

M = x(1)*M_AA + x(2)*M_AB + x(3)*M_BA + x(4)*M_BB;
IM = eye(5)/M;
%IM = pinv(M);
y = log(IM(3,3));

