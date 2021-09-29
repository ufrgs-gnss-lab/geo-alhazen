function psi = z2psi (zs, phi, k1, k2)
  
% Evaluation for the correct root

% Roots of quartic polynomial ([alpha, beta, gamma, conj(beta), conj(alpha)])
psis = (1/2)*acos(real(zs));  % bottom p.473

% Condition: phi + 2*psis = acos(k1*cos(psis)) +acos(k2*cos(psis)) % eq.2, p.473
% Left hand condition
eq2_lhs = phi +2*psis;
% Right hand condition
eq2_rhs = acos(k1*cos(psis)) +acos(k2*cos(psis)); % eq.2, p.473
% Equation 2 equality test
eq2_neq = eq2_lhs - eq2_rhs;

% Column of mininum eq2_neq (Equality test)
eq2_col = argmin(abs(eq2_neq), [], 2);

[n, m] = size(zs);

% [n,m]Size of array, (1:n)Rows of array, eq2_col column of array
% Return the indices
eq2_ind = sub2ind([n m], (1:n)', eq2_col);
% Get the minimum error of equalities test
eq2_err = eq2_neq(eq2_ind);
% Tolerance
tol = sqrt(eps());
% Error should be smallest than tolerance
eq2_cond = (abs(eq2_err) < tol);
% "the unique value psi of psi* that satisfies [eq.(2)] is the grazing angle" bottom p.473
% Psi receive the real part of correct root of quartic polynomial
psi = real(psis(eq2_ind));

end
