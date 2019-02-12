% project name  : Ising Model 
% Author        : Hosein Ghorban Fekr - Hoseinqf@yahoo.com
% Creation Date : 2009/04/12
% Description   : This program simulates Ising Model in 2D
% --------------------------------------------------------
clc
clf
clear
global S L;
% --------------
L   = 50;   % 2D system L*L
Ns  = 1000;  % no. of ensembles
% -------------------------------------------------
N   = L^2;
M   = 0;
E   = 0;
Stp = 0;
S   = zeros(L,L);
for n=1:N  S(n) = rnd; M = M + S(n); end;
for i=1:L-1
    for j=1:L-1
        E = E + 2 * S(i,j) * (S(i+1,j) + S(i,j+1));
    end
end
% -------------------------------------------------
t = 0;
for T=5:-0.1:0
    k   = 0;
    Bet = 1/(T+eps);
    fprintf('T = %1.2f ...',T);
    a   = [exp(8*Bet) exp(4*Bet) 1 exp(-4*Bet) exp(-8*Bet)];
    for n=1:Ns*N
        i   = fix(rand*L+1);
        j   = fix(rand*L+1);
        tmp = DelE(i,j);      % energy diffrerence
        if tmp<0
            S(i,j) = -1 * S(i,j);
            E      = E + tmp;
            M      = M + 2 * S(i,j);
        elseif a(tmp/4+3)>rand
            S(i,j) = -1 * S(i,j);
            E      = E + tmp;
            M      = M + 2 * S(i,j);
        end
        % -----------------------
        if mod(n,N)==0
            k = k + 1;
            En(k) = E;
            Mn(k) = M;
%             if Stp==0 && mod(n,5*N)==0
%                  pltS;         % animate system
%                 if abs(sum(sum(S(:,1:L))))==N  Stp = 1; end
%             end
        end
    end
    % ------------------
    t = t + 1;
    N   = length(En);
    if t==1
        data = [T mean(En) mean(Mn) Bet^2*var(En) Bet*var(Mn)];
    else
        data = [data;T mean(En) mean(Mn) Bet^2*var(En) Bet*var(Mn)];
    end
    clc
    clear En Mn;
end
save data.txt data -ascii
subplot(221)
plot(data(:,1),data(:,2),'.')
title('Energy - T')
subplot(222)
plot(data(:,1),abs(data(:,3)),'.')
title('M - T')
subplot(223)
plot(data(:,1),data(:,4),'.')
title('C - T')
subplot(224)
plot(data(:,1),data(:,5),'.')
title('Kapa - T')