function J = PI(K, A, B, C, Q, R)

A_c = A - B*K*C;

max_real_eig = max(real(eig(A_c)));

if max_real_eig < 0
    P = lyap(A_c.', Q+C.'*K.'*R*K*C);
    J = 0.5*trace(P);

else
    J = 1e10+1e20*max_real_eig;

end
