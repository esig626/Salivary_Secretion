
clear

% close all

clc

format longg

float('double')
e_i = 1676.42325017681;
e_l = 5616.28349757687;


IC = [0.982954608922101
        0.0156056873998078
        0.0698613121949313
        0.0275833616327992
       0.00570251673610959
        0.0212577961855938
      0.000841285060578805
        0.0797272597523379
                         0
         0.503415009836623
         0.779722071043159
        0.0595238376705087];

IC = [0.885986458183828
        0.0234646251685931
        0.0662964483773376
        0.0209743524206292
       0.00509993352666838
        0.0227747944603051
      0.000931251810034654
         0.103516197255918
        0.0142355306859098
        0.0303860755935395
         0.666203550181955
         3];
     
IC = double(IC)

V_PLC = 0.0013;

f = @(t,x) ND(t,x,V_PLC);
g = @(t,x) ND_OPEN(t,x,V_PLC);

t_end = 600;

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[t,U_init] = ode15s(f,[0,t_end],IC,options);
[t,U] = ode15s(f,[0,300],U_init(end,:),options);
[t_OPEN,U_OPEN] = ode15s(g,[0,60],U_init(end,:),options);

v   = U(:,1);            
xi  = U(:,2);           
yi  = U(:,3);            
zi  = U(:,4);           
bi  = U(:,5);                     
xl  = U(:,6);           
yl  = U(:,7);
ca  = U(:,8);
ip  = U(:,9);
h   = U(:,10);
g   = U(:,11);
ct  = 3;
%hi  = U(:,13);


            
Na  = e_i *U(:,2);           
K  = e_i *U(:,3);            
Cl  = e_i *U(:,4);           
HCO3  = e_i *U(:,5);                     
Nal  = e_l *U(:,6);           
Kl  = e_l *U(:,7);

XV = 46.2653795533866;
Hy = - (xi + yi - zi - bi - XV./(v*904))*e_i;



%xv = (v*904).*((10^(-6)*hi)+xi+yi-bi-zi);


hy = - (xi + yi - zi - bi - XV./(v*904)).*10^6;


ylab = {'Volume','Na','K','Cl','HCO3','Nal','Kl','Ca','IP_3','h','g','c_t'};

scale = [904 e_i e_i e_i e_i e_l e_l 1 1 1 1 1];
param.Cle   =102.6;          %mM
param.Nae   =140.2;          %mM
param.Ke    =5.3;             %mM
param.Hel    =3.8904514499428e-05;     %mM
param.HCO3e =40;
param.La = 4.43354480372293;
param.Lb = 1.94580852061638;
param.Lt = 0.0593663800587253;
param.Ul = 51.7310374376987;   
param.CO20 = 6.60012219974207;
Qa = param.La * ( 2 * ( Nal + Kl - Na - K - Hy ) - param.CO20 + param.Ul ); 

Qb = param.Lb * ( 2 * ( Na + K + Hy ) + param.CO20 - ...
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );

Qt = param.Lt * ( 2 * ( Nal + Kl ) + param.Ul - ....
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); 

Qtot=(Qa+Qt);
h= figure(3)
plot(t,Qtot,'LineWidth',3.5,'Color','k')
% hold on
% plot(t_OPEN,Qtot_OPEN,'LineWidth',3.5,'Color','k')
xlabel('Time (s)')
ylabel('Fluid flow \mum^3.s^-1')
box off
ax=gca;
ax.LineWidth = 10;
ax.FontSize = 40;
xlim([0,62])
xticklabels({'0','20','40', '60'})
xticks([0 20 40 60])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
hold on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,['FF_ODE.pdf'],'-dpdf','-r0')
figure(1)
for i = 1:12
    subplot(6,2,i)
    plot(t,scale(i)*U(:,i),'LineWidth',2);
    hold on
    ax = gca;
    ax.LineWidth = 2.1;
    ax.FontSize = 20;
    xlabel('Time (sec)')
    ylabel(ylab{i})
    box off
end
% 
% subplot(2,2,1)
% plot(t,ca,'LineWidth',3);
% hold on
% ax = gca;
% ax.LineWidth = 2.1;
% ax.FontSize = 20;
% xlabel('t')
% ylabel('[Ca^{2+}] (\muM)')
% box off
% subplot(2,2,2)
% plot(t,zi*e_i,'LineWidth',2)
% hold on
% ax=gca;
% ax.LineWidth = 3; 
% ax.FontSize = 20;
% xlabel('Time (sec)')
% ylabel('[Cl^-] (mM)')
% box off
% subplot(2,2,3)
% plot(t,yi*e_i,'LineWidth',2)
% hold on
% ax=gca;
% ax.LineWidth = 3; 
% ax.FontSize = 20;
% xlabel('Time (sec)')
% ylabel('[K^+] (mM)')
% box off
% subplot(2,2,4)
% plot(t,hy,'LineWidth',2)
% hold on
% ax=gca;
% ax.LineWidth = 3; 
% ax.FontSize = 20;
% xlabel('Time (sec)')
% ylabel('Hydrogen')
% box off

h = figure(5)
plot(t,U(:,8),'LineWidth',3.5,'Color','k')
% hold on
% plot(t_OPEN,U_OPEN(:,8),'LineWidth',3.5)
ax=gca;
ax.LineWidth = 10;
ax.FontSize = 100;
xlabel('Time (sec)')
ylabel('[Ca^{2+}] (\muM)')
box off
xlim([0,302])
xticklabels({'0','100','200', '300'})
xticks([0 100 200 300])
ylim([0.0788, .45])
% yticks([0.0789 0.2 0.27])
% yticklabels({'70','200','270'})
% text(65,0.395,...
%         '\underline{\qquad \qquad Agonist \qquad \qquad}', ...
%         'FontSize',95,'Interpreter', 'latex','FontName','Helvetica')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
hold on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,['Ca_ODE.pdf'],'-dpdf','-r0')
