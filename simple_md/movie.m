% Hossein Ghorban Fekr - 87206975 - hoseinqf@yahoo.com
% molecular Dynamic simulating program
% This program read data from main program which wrote in C++. 
clc
clf
clear
NL = [100 100];
R = load('RV2.dat');
tic
tmp=toc;
for j=0:length(R)/NL(1)-1
        while 1
            if toc-tmp>0.001
                plot(R(j*NL(1)+1:(j+1)*NL(1),1),R(j*NL(1)+1:(j+1)*NL(1),2),'o','markersize',3,'markerfacecolor','b')
                axis([0 NL(2) 0 NL(2)]);
                drawnow;
                tmp = toc;
                break;
            end
        end
end


%  d= load('UKETPn.dat');
%  hold on
%  plot(1:length(d),d(:,6),'o','markersize',1,'markerfacecolor','b')
%  plot(1:length(d),d(:,2),'o','markersize',1,'markerfacecolor','b')
%  plot(1:length(d),d(:,3),'o','markersize',1,'markerfacecolor','b')
%  drawnow
%  waitfor(1)