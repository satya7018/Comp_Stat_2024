function LL = likelihood_Frank(theta,ns)
%theta = [1 2 0.5 0.4 2];
x = [0.25 0.25 0.25 0.25];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n11_AA=ns(1);
n01_AA=ns(2);
n10_AA=ns(3);
n00_AA=ns(4);

n11_AB=ns(5);
n01_AB=ns(6);
n10_AB=ns(7);
n00_AB=ns(8);

n11_BA=ns(9);
n01_BA=ns(10);
n10_BA=ns(11);
n00_BA=ns(12);

n11_BB=ns(13);
n01_BB=ns(14);
n10_BB=ns(15);
n00_BB=ns(16);

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

C_AA = -(1/theta(5))*log(1+(((exp(-theta(5)*pi1_AA)-1)*(exp(-theta(5)*pi2_AA)-1))/(exp(-theta(5))-1)));
C_AB = -(1/theta(5))*log(1+(((exp(-theta(5)*pi1_AB)-1)*(exp(-theta(5)*pi2_AB)-1))/(exp(-theta(5))-1)));
C_BA = -(1/theta(5))*log(1+(((exp(-theta(5)*pi1_BA)-1)*(exp(-theta(5)*pi2_BA)-1))/(exp(-theta(5))-1)));
C_BB = -(1/theta(5))*log(1+(((exp(-theta(5)*pi1_BB)-1)*(exp(-theta(5)*pi2_BB)-1))/(exp(-theta(5))-1)));

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

%l_AA= y1*y2*log(p11_AA)+y1*(1-y2)*log(p10_AA)+(1-y1)*y2*log(p01_AA)+(1-y1)*(1-y2)*log((1-p11_AA-p10_AA-p01_AA));
%l_AB= y1*y2*log(p11_AB)+y1*(1-y2)*log(p10_AB)+(1-y1)*y2*log(p01_AB)+(1-y1)*(1-y2)*log((1-p11_AB-p10_AB-p01_AB));
%l_BA= y1*y2*log(p11_BA)+y1*(1-y2)*log(p10_BA)+(1-y1)*y2*log(p01_BA)+(1-y1)*(1-y2)*log((1-p11_BA-p10_BA-p01_BA));
%l_BB= y1*y2*log(p11_BB)+y1*(1-y2)*log(p10_BB)+(1-y1)*y2*log(p01_BB)+(1-y1)*(1-y2)*log((1-p11_BB-p10_BB-p01_BB));

l_AA= n11_AA*log(p11_AA)+n10_AA*log(p10_AA)+n01_AA*log(p01_AA)+n00_AA*log((1-p11_AA-p10_AA-p01_AA));
l_AB= n11_AB*log(p11_AB)+n10_AB*log(p10_AB)+n01_AB*log(p01_AB)+n00_AB*log((1-p11_AB-p10_AB-p01_AB));
l_BA= n11_BA*log(p11_BA)+n10_BA*log(p10_BA)+n01_BA*log(p01_BA)+n00_BA*log((1-p11_BA-p10_BA-p01_BA));
l_BB= n11_BB*log(p11_BB)+n10_BB*log(p10_BB)+n01_BB*log(p01_BB)+n00_BB*log((1-p11_BB-p10_BB-p01_BB));

LL = -(x(1)*l_AA + x(2)*l_AB + x(3)*l_BA + x(4)*l_BB);
end