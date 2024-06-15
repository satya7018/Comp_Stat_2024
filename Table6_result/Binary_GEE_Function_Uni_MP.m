function y = Binary_GEE_Function_Uni_MP(ps,X,rho,n)
%X = lhsdesign(100,5);
%s = [];
% FC R
% lb = [-3.157 -1.986 -4.323 -7.298]; 
%  ub = [3.772 1.978 3.701 6.557];

% GC R
% lb = [-3.154 -1.985 -4.224 -7.102]; 
%  ub = [3.754 1.974 3.70 6.554];

%  % CC R
% lb = [-3.154 -1.985 -4.224 -7.102]; 
%  ub = [3.754 1.974 3.700 6.554];

 % NC R
% lb = [-3.160 -1.987 -4.303 -7.262]; 
%  ub = [3.760 1.974 3.708 6.574];

  % FC H3
lb = [-4 0 0 0]; 
 ub = [2 3 2 2];

%  lb = [-3.154 -1.985 -4.224 -7.102]; % lower bounds for A,V,h,l and b
%  ub = [3.754 1.974 3.700 6.554];
% lb = [0 0 0 0]; % lower bounds for A,V,h,l and b
% ub = [0.9324 0.56 1.3873 1.1553];

T = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));
%  ps = [1 0 0 0];
%  rho = 0.0;
for i = 1:n
    s(i) = Binary_GEE_Function(ps,T(i,1),T(i,2),T(i,3),T(i,4),rho);
end
% for i = 1:100
%     s2(i) = p2t2bal(T(i,1),[T(i,2) T(i,3)],T(i,4),T(i,5),rho);
% end
y = ((mean(s)));
%y(2) = ((mean(s2)))