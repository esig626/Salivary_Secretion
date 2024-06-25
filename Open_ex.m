
run Int_par

close all

V_PLC = linspace(0,0.0013,100);
ALPHA_0 = [0.0027, 0.0027 * 2];
colour = ['b','r'];

for k = 1:length(ALPHA_0)
for j = 1:length(V_PLC)

g = @(t,x) ND_OPEN(t,x,V_PLC(j),ALPHA_0(k));

[t,U] = ode15s(g,[0,t_end],IC,options);

v   = U(:,1);            
xi  = U(:,2);           
yi  = U(:,3);            
zi  = U(:,4);           
bi  = U(:,5);                     
XV = 46.2653795533866;
Hy = - (xi + yi - zi - bi - XV./(v*904))*e_i;


for i = 1:11 
   U(:,i) = scale(i)*U(:,i);
end

v   = U(:,1);            
Na  = U(:,2);           
K  = U(:,3);            
Cl  = U(:,4);           
HCO3  = U(:,5);                     
Nal  = U(:,6);           
Kl  = U(:,7);
ca  = U(:,8);
ip  = U(:,9);
h   = U(:,10);
g   = U(:,11);
ct  = U(:,12);


% figure(1)
% plot(t,ca,'LineWidth',2)
% ax = gca;
% ax.LineWidth = 2;
% ax.FontSize = 15;
% box off
% ylabel('Cytosolic [Ca^{2+}]')
% xlabel('Time')
% hold on


Qa = param.La * ( 2 * ( Nal + Kl - Na - K - Hy ) - param.CO20 + param.Ul ); 

Qb = param.Lb * ( 2 * ( Na + K + Hy ) + param.CO20 - ...
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );

Qt = param.Lt * ( 2 * ( Nal + Kl ) + param.Ul - ....
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); 

Qtot=(Qa+Qt);

M_Q(1,j) = mean(Qtot);

% figure(2)
% plot(t,Qtot,'LineWidth',2)
% hold on
% ax = gca;
% ax.LineWidth = 2;
% ax.FontSize = 15;
% box off
% ylabel('Flow Rate')
% xlabel('Time')
% legend(num2str(M_Q(1,j)))
% legend boxoff


figure(1)
plot(V_PLC(1,j),M_Q(1,j),'.','Color',colour(k),'MarkerSize',15)
hold on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 15;
box off
ylabel('Flow Rate')
xlabel('V_{PLC}')


end
end