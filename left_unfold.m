function LU = left_unfold(U)
[l,c,r] = size(U);
LU = reshape(U,[l*c, r]);
end