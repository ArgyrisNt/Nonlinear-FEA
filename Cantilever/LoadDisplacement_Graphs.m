%% Plot graphs for the initial structure and for the modified one

clear all;
clc;

% Read data for loads,displacements for the initial structure
fileID1 = fopen('ForcesDisplacementsData(1).txt','rt');
data1 = textscan(fileID1, '%d %f', 'HeaderLines',2);
fclose(fileID1);
P1 = data1{1,1};
u1 = data1{1,2};

% Read data for loads,displacements for the modified structure
fileID2 = fopen('ForcesDisplacementsData(2).txt','rt');
data2 = textscan(fileID2, '%d %f', 'HeaderLines',2);
fclose(fileID2);
P2 = data2{1,1};
u2 = data2{1,2};

% Plot the two load-displacement curves
tiledlayout(2,1)
nexttile
plot(u1,P1)
grid on
title('Plot 1')

nexttile
plot(u2,P2)
grid on
title('Plot 2')