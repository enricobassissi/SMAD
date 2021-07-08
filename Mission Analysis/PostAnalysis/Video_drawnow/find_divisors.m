function D = find_divisors(n)

K=1:n
D = K(rem(n,K)==0)

end