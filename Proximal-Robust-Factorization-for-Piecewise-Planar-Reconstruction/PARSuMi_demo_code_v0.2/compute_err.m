function err=compute_err(N,M,C,mask)
Wr=N*C;
err=0.5*norm(M-Wr(mask))^2;
end