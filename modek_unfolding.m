function mat = modek_unfolding(ten, n)

    N = ndims(ten);
    order = [n n+1:N 1:n-1];
    perm_ten = permute(ten, order);
    ten_size = size(perm_ten);
    rows     = ten_size(1);
    cols     = prod(ten_size(2:end));
    mat      = reshape(perm_ten, [rows, cols]);
    % mat      = mat.data; 
    
    
end