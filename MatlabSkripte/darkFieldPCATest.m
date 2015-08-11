load hald;

data = ingredients(:,1:3);

data = bsxfun(@minus,data,mean(data));

data

figure(1);
plot3(data(:,1),data(:,2),data(:,3),'.');
title('Three dimensional plot');

[pc,score,latent,tsquare] = princomp(data);

pc

figure(2);
biplot(pc(:,1:2),'Scores',score(:,1:2),'VarLabels',...
		{'X1' 'X2' 'X3'})