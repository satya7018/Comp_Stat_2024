function y = BCD_product_Uniform(xx,lb,ub)
n=2000;
%mu=[0.5,0.5,0.5,0.5,0.5];
%sigma= [1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,1];
X = lhsdesign(n,4);
%X1 = lhsnorm(mu,sigma,n);
%s = [];
% lb = [1, -1, 1, -1]; % lower bounds for A,V,h,l and b
% ub = [2, 1, 2, 1];
T = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));

for i = 1:n
    %s(i) = Normal_Inf_Matrix([T(i,1),T(i,2),T(i,3),T(i,4),T(i,5)],xx);
     s(i) = Product_Inf_Matrix([T(i,1),T(i,2),T(i,3),T(i,4)],xx);
end

y = ((mean(s)));

end