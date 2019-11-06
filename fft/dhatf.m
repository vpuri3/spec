function Dhf=dhf_mat(Mi);
%
%  Fourier derivative matrix for sine/cosine bases
%

M=Mi; if mod(M,2)==0; M=M+1; end; M1=M-1;  %% M1 is even

Dhf=zeros(M+1,M+1);

K=diag(0:M1/2);
S=[0 1;-1 0]; S=sparse(S);
Dhf=kron(K,S);            %% Dhf is sparse because S is sparse
Dhf=Dhf(2:end,2:end);     %% Get rid of (0,0) mode
Dhf=Dhf(1:Mi,1:Mi);       %% Retain only Mi entries



