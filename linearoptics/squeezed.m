function state = squeezed(z,N)
%%% Function that calculates the number state coefficients of a squeezed state
%%% Inputs
% z : complex, squeezing coefficient
% N : the maximum of total photons (for Fock space truncation)
%%% Outputs
% state : an array storing the coefficients of |0>, |1>, |2>, ..., |N>

r = abs(z);
t = tanh(r);
state = zeros(1,N+1);
a = -0.5*(z/r)*t;
for itN = (1:floor(N/2)+1)
    state(2*itN-1) = a^(itN-1)*sqrt(factorial(2*(itN-1)))/factorial(itN-1);
end
state = state/sqrt(cosh(r));

end