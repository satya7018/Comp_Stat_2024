

rngd1 = randi(1000)

rng(rngd1)

%rng(613)

% GEE based optimal design
lb1 = [0 0 0 0];
ub1 = [1 1 1 1];
x1 = [0.218 0.344 0.168 0.270];
nvars = 4;
x11 = [];
fval1 = [];
n = 100;
X = lhsdesign(n,4);
% Uniform prior range (a,b)
a = 0;
b = 0.5;

options = optimoptions('fmincon','Algorithm','interior-point','TolFun',1e-8,'TolCon',1e-8,'MaxFunEvals',10000);
[x11, fval1] = fmincon(@(ps)Binary_GEE_Function_Uni_MP_Uniform(ps,X,n,a,b),x1,[],[],[1 1 1 1],1,lb1,ub1,[],options)
%x11 will be GEE based design 

% FC H3
 ps11 = [0.235 0.308 0.171 0.286];

%  % GC H3

ps12 = [0.218 0.344 0.168 0.270];

% % CC H3
  ps13 = [0.280 0.332 0.141 0.247];

% NC H3
  ps14 = [0.233 0.321 0.178 0.268];

s11=[];
s21=[];
es11=[];

s12=[];
s22=[];
es12=[];

s13=[];
s23=[];
es13=[];

s14=[];
s24=[];
es14=[];


for i = 1:1000
    %rho = betarnd(a,b);
        rho = a + (b-a)*rand(1);
s11(i) = Binary_GEE_Function_Uni_MP(ps11,X,rho,n);
s21(i) = Binary_GEE_Function_Uni_MP(x11,X,rho,n);
es11(i) = s11(i)/s21(i);

s12(i) = Binary_GEE_Function_Uni_MP(ps12,X,rho,n);
s22(i) = Binary_GEE_Function_Uni_MP(x11,X,rho,n);
es12(i) = s12(i)/s22(i);

s13(i) = Binary_GEE_Function_Uni_MP(ps13,X,rho,n);
s23(i) = Binary_GEE_Function_Uni_MP(x11,X,rho,n);
es13(i) = s13(i)/s23(i);

s14(i) = Binary_GEE_Function_Uni_MP(ps14,X,rho,n);
s24(i) = Binary_GEE_Function_Uni_MP(x11,X,rho,n);
es14(i) = s14(i)/s24(i);
end
[mean(es11) mean(es12) mean(es13) mean(es14)]


 

