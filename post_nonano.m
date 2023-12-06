clc; clear; close all;

load dados.mat

Tm = u;
l = length(Tm);
T_max = max(max(Tm));
[Tmax_idxx, Tmax_idxy] = find(Tm==T_max);

T_x = Tm(250,:);

[~,TT] = find(T_x>41,1);

[~, ridx] = min(abs(x-0.055));
[~, ridy] = min(abs(y-0.055));

x_163 = x(163);
x_184 = x(184);
xHealthy = round(x_184 - x_163,3);

T184 = Tm(250,163)
r = 2;
R = 2.6;
areaSaudavel = pi*(R^2-r^2);
Ar = 12.57;
AR = 21.24;
rDiff = 100*(Ar-AR)/AR

% hfig = figure;
% plot(x,T_x,'LineWidth',1.5)
% axis square;
% xlabel('Comprimento (m)')
% ylabel('Largura (m)')
% xline(5.5/100,'k--','Limite do tumor','LineWidth',1)
% xline(9.5/100,'k--','Limite do tumor','LineWidth',1)
% yline(41,'r--','Necrose celular')
% fname2 = 'linegraph_nonano';

% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',14) % adjust fontsize to your document

% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

% print(hfig,fname2,'-dpdf','-vector','-bestfit')