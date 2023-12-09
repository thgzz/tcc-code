clc; clear; close all;

load dados_nanor11.mat

% vetor tempo 300s
tp = linspace(0,T,500);

Tg = u2;
T_max = max(max(Tg))
[Tmax_idxx, Tmax_idxy] = find(Tg==T_max);

Txx = Tg(250,:);

[~,TT] = find(Txx>41,1)
rr = 2;
x182 = x(182)
x184 = x(184)
x250 = round(x(250),3);
xfinal = (x184-x182)*100;
xx = (x250*100)+(rr+xfinal);
ff = round(xx*100,2);

[~, ridx] = min(abs(x-0.055));

T184 = Tg(250,184)
T182 = Tg(250,182)

r_187 = round(x(187),3)
