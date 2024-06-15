function y = BCD_Clayton_Normal(xx, mu, lb, ub, sgm)
n=1000;
%mu=[0.5,0.5,0.5,0.5,0.5];
%sigma= [1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
X = lhsnorm(mu,sgm*eye(4),n);
r = [];
for i = 1:n
    r(i) = unifrnd(lb,ub);
end

T = [X r'];
% ps = [0.25 0.25 0.25 0.25];
% rho = 0.2;
for i = 1:n
    s(i) = Clayton_Inf_Matrix([T(i,1),T(i,2),T(i,3),T(i,4),T(i,5)],xx);
end
% for i = 1:100
%     s2(i) = p2t2bal(T(i,1),[T(i,2) T(i,3)],T(i,4),T(i,5),rho);
% end
y = ((mean(s)));
%y(2) = ((mean(s2)))
end