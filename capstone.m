% ms = [1 1];
% len = 2;
 ms = [.058 .548 .201 .101];
 len = [.1651 .2032 .1778 .1778];
 mass = diag(ms);
 length = diag(len);
% n = 4;         % number of masses
   n = 4;
 g = 9.8;       % gravitational constant
 %g = 1;
 pi = 3.1415;
% set up system matrices

% Fix the off diagonal entries in A !!! (Kenny Davis)

%A =  (g*m)/(2*l)*(diag(1:2:2*n-1)-diag(1:n-1,1)-diag(1:n-1,-1)); % correct
K = (g*mass)./(2.*len)*(diag(1:2:2*n-1)-diag(1:n-1,1)-diag(1:n-1,-1));
xdouble = -1*(inv(mass)\K).*len; %invert mass?

 % compute modes U and eigenvalues diag(D)
 [U, W2] = eig(mass\K); 

% sort modes so the ones with lowest frequency appear first
[junk,indx] = sort(diag(W2));
U = U(:,indx);
W2 = W2(indx,indx);
U = U*diag(1./U(n,:));
f2 = W2/(4*pi^2);
f = sqrt(f2);