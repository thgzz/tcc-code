clc; clear; close all;


load dados_nano.mat
load dados_nanor11.mat

% vetor tempo 300s
tp = linspace(0,T,500);

Tm = u;
Tg = u2;
l = length(Tm);
% [T_max, idx] = max(Tm,[],"all")
T_max = max(max(Tm))
[Tmax_idxx, Tmax_idxy] = find(Tm==T_max);

T_x = Tm(250,:);
Txx = Tg(250,:);

[~,TT] = find(T_x>41,1);

[~, ridx] = min(abs(x-0.055));

T184 = Tm(250,184)

r_187 = round(x(187),3)
R = 2;
r = 1.9;
anlA = pi*(R^2-r^2)

AR = 12.57;
Ar = 11.34;

rDiff = 100*(AR-Ar)/Ar;


hfig = figure;
plot(x,T_x,'LineWidth',1.5)
hold on
plot(x,Txx,'LineWidth',1.5)
hold off
legend('$r_0 = 10 \;mm$','$r_0 = 11 \;mm$')
% axis square;
xlim([0 0.15]);
ylim([30 110])
xticks([0:0.05:0.15]);
xtickformat('%.2f');
xlabel('Comprimento (m)')
ylabel('Temperatura (°C)')
xline(5.5/100,'k--','Limite do tumor','LineWidth',1,'HandleVisibility','off')
xline(9.5/100,'k--','Limite do tumor','LineWidth',1,'HandleVisibility','off')
yline(41,'r--','Necrose celular','HandleVisibility','off')
fname2 = 'linegraph_nano';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document

set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

print(hfig,fname2,'-dpdf','-vector','-bestfit')