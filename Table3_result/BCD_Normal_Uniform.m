function y = BCD_Normal_Uniform(xx,lb,ub)
n=1000;
X = lhsdesign(n,5);
T = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));
s =[];
for i = 1:n
    s(i) = Normal_Inf_Matrix([T(i,1),T(i,2),T(i,3),T(i,4),T(i,5)],xx);
end

y = mean(s);

end