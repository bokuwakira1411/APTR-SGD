function Z = Initialize(T,r)
I = size(T);
d = length(I);
Z = {};
Z{1} = randn(r(d), I(1), r(1));
for i=2:1:d
    Z{i} = randn(r(i-1),I(i),r(i));
end
end
