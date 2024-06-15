function y = BCD_Frank_Normal(xx, mu, lb, ub, sgm)
n=1000;
X = lhsnorm(mu,sgm*eye(4),n);
r = [];
for i = 1:n
    r(i) = unifrnd(lb, ub);
end

T = [X r'];

for i = 1:n
    s(i) = Frank_Inf_Matrix([T(i,1),T(i,2),T(i,3),T(i,4),T(i,5)],xx);
end
y = mean(s);

end