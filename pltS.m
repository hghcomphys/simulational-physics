% pltS
global S L;
[row,col] = find(S==1);
% subplot(121)
plot(row,col,'.','Markersize',20,'color','g');
axis([1 L 1 L]);
axis square
[row,col] = find(S==-1);
hold on;
plot(row,col,'.','Markersize',20,'color','b');
hold off
drawnow