function RU = right_unfold(U)
[l,c,r] = size(U);
RU = reshape(U, [l, c*r]);
end