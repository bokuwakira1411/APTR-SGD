function T = modek_classical_unfolding(U,k)
I = size(U);
d = length(I);
set = [k 1:k-1 k+1:d];
T = permute(U,set);
T = reshape(T,[I(k), prod(I)/I(k)]);
end
