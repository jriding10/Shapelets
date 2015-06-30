%% Create hermite polynomial matrix

n_max = 100;
final = zeros(n_max+1);
for n = 0:n_max
    builder = hermite_rec(n);
    for k = 1:(n+1)
        final(n+1, k) = builder(k);
    end
end

save('hermite_n100_coeffs.txt', 'final', '-ascii', '-double');

