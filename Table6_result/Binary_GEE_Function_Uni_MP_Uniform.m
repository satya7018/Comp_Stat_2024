function y = Binary_GEE_Function_Uni_MP_Uniform(ps,T,n,a,b)
for i = 1:100
    rho = a + (b-a)*rand(1);
s(i) = Binary_GEE_Function_Uni_MP(ps,T,rho,n);
end
y = mean(s);