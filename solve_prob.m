function [x1, zi] = solve_prob(lambda, Q, r, K, xlow, xup,n )
cvx_begin 
cvx_solver mosek
cvx_solver_settings('MSK_DPAR_MIO_MAX_TIME',20)
variables x1(n) 
variable zi(n) binary
%                                  Object function
minimize ( (lambda *(x1'*Q*x1)) -(1-lambda)*(sum(r'*x1)))
%                                  Constrait(s)
subject to
sum(x1)<=1; %Proportion of the total funds invested in asset is equal to 1
sum(zi)==K; % Cardinality constraint that limits the number of assets
x1 >= xlow*zi;
x1 <= xup*zi;
x1 >= 0;
cvx_end
end

