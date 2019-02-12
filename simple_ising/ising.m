% Name  : Ising Model 
% Author        : Hosein Ghorbanfekr - hgh.comphys@gmailc.com
% Creation Date : 2009/04/12
% Description   : This program simulates Ising Model in 2D
%                 at given temperature in order to study its
%				  the critical behavior.
% --------------------------------------------------------
clc; clf; clear;
global S L;
% --------------
L   = 50;   % 2D system L*L
T   = 1; % Temprature
Ns  = 300;  % no. of ensembles
Nc  = 75;   % fit point for finding To.
% --------------
N   = L^2;
Bet = 1/T;
S   = zeros(L,L);
Stp = 0;
a   = [exp(8*Bet) exp(4*Bet) 1 exp(-4*Bet) exp(-8*Bet)];
for n=1:N  S(n) = rnd; end; X(:,:,1) = S;
% ------------------------------------------------------
k = 1;
for n=1:Ns*N
    i   = fix(rand*L+1);
    j   = fix(rand*L+1);
    tmp = DelE(i,j);      % energy diffrerence
    if tmp<0
        S(i,j) = -1 * S(i,j);
    elseif a(tmp/4+3)>rand
        S(i,j) = -1 * S(i,j);
    end
    % -----------------------
    if mod(n,N)==0
        k = k + 1;
        X(:,:,k) = S; % saving ensembles in X.mat
        if Stp==0
            pltS;         % animate system
            if abs(sum(sum(S(:,1:L))))==N  Stp = 1; end
        end
    end
end
save X;
% ////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\
N = size(X,3);
E = zeros(N,1);
M = zeros(N,1);
for n=1:N
    % ------------------
    for i=1:L-1
        for j=1:L-1
            E(n) = E(n) - X(i,j,n)*(X(i+1,j,n)+X(i,j+1,n)); % calculate E
        end
    end
    % ------------------
    M(n) = sum(sum(X(:,1:L,n))); % calculate M
end
% ---------------------------------
C   = zeros(N,1);
tmp = var(E);
for j=0:N-1
    a = [ 0 0 0 ];
    for i=1:(N-j)
        a(1) = a(1) + E(i)*E(j+i);
        a(2) = a(2) + E(i);
        a(3) = a(3) + E(i+j);
    end
    a      = a/(N-j);
    C(j+1) = (a(1)-a(2)*a(3))/tmp; % calculating Cj
end
subplot(122)
plot(0:N-1,C,'.','color','r') % plot C - j graph
ylabel('C(j)');
xlabel('j');
f = fittype('a*exp(-x/k)+b');
F = fit((0:Nc-1)',C(1:Nc),f); clc %fit graph for finding kisi
T = fix(F.k)+1;
% --------------------- print data ----------------------
fprintf('To = %d\n',T);
fprintf('E = %f , ErrE = %f\n',mean(E),std(E)/sqrt(N/T));
fprintf('M = %f , ErrM = %f\n',mean(M),std(M)/sqrt(N/T));
fprintf('C = %e\n',Bet^2*var(E));
fprintf('K = %e\n\n',Bet*var(M));
