clc; clear; close all;


load dados_nano.mat

% vetor tempo 300s
tp = linspace(0,300,250);

Tm = u;
l = length(Tm);
[T_max,idx] = max(Tm,[],"all")
M = max(Tm,[],2);

T_x = Tm(250,:);

T_idx = find(T_x>37,1)

plot(x,T_x)
xline(5.5/100,'k--')
xline(9.5/100,'k--')
yline(41,'r--')