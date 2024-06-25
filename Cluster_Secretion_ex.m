clear

close all

clc

format longg

load('PSec.mat')


c7 = Cell(7,[6,5],psec);
c6 = Cell(6,[7,5,4],psec);
c5 = Cell(5,[7,6,3,2],psec);
c4 = Cell(4,[6,2,1],psec);
c3 = Cell(3,[5,3,2],psec);
c2 = Cell(2,[5,4,3,1],psec);
c1 = Cell(1,[4,2,1],psec);

%%

h = @(t,x) Cluster_Secretion(t,x,c7,c6,c5,c4,c3,c2,c1);

IC = [915.723654477861,24.1909570454928,120.731778907326,...
     51.6694753355558,...
     10.0972360642965,118.730344586928,4.53763537277269,...
     169.699581870451,...
     118.733409170121,4.54483829049721,114.064438951154,...
     855.541686626341,24.3241213203340,120.826663018036,...
     51.4248489435880,10.0816305795372,118.696293282314,...
     4.52431428484669,112.050734621348,118.683043979139,...
     4.51622904664592,144.441847432793,374.955401545744,...
     30.7171320033499,114.594015477197,43.6370739575270,...
     9.47800982229926,118.733982004301,4.52065065348946,...
     78.8202193782775,118.716737215782,4.51491865511491,...
     116.539464845002,873.092918436770,24.6732432125689,...
     120.499448882435,50.9140501586920,10.0414361456220,...
     118.742018290914,4.51267747170031,120.813399852893,...
     118.737692579792,4.51617776948688,123.253870349272,...
     652.902993316093,26.3375576398869,118.964241413937,...
     48.6165858420678,9.86305992011740,118.733418708820,...
     4.51264074240377,82.2720472345755,798.450950914157,...
     24.5637976487375,120.570918055794,51.0796131198430,...
     10.0539261572363,118.707242076041,4.51485139054405,...
     155.725041554827,118.700910598084,4.51032530510424,...
     124.064497069754,791.901307068352,25.2029146936759,...
     119.952101326014,50.1699727610863,9.98236956014850,...
     118.717385697454,4.51687419398050,116.806864108262,...
     0.0703634335784682,0,0.625605708678142,0.819641890101789,...
     0.0702891893343451,0,0.626594224643312,0.819953808976032,...
     0.0700955124233116,0,0.629173837092982,0.820767058954173,...
     0.0703562624009537,0,0.625701312714952,0.819672021981437,...
     0.0701528556794530,0,0.628409476631628,0.820526340774401,...
     0.0705241664847193,0,0.623465658084026,0.818966296960795,...
     0.0700640546820776,0,0.629591904352207,0.820899091033011];

t_end = 500;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-16);

[t0,U] = ode15s(h,[0,t_end],IC,options);

[t,x] = ode15s(h,[0,t_end],U(end,:),options);


clear options h t_end IC

%%
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%
w_7   = x(:,1); xi_7  = x(:,2); yi_7  = x(:,3); zi_7  = x(:,4); 
bi_7  = x(:,5); xl_76 = x(:,6); yl_76 = x(:,7); zl_76 = x(:,8); 
xl_75 = x(:,9); yl_75 = x(:,10); zl_75 = x(:,11); w_6 = x(:,12); 
xi_6  = x(:,13); yi_6  = x(:,14); zi_6  = x(:,15); bi_6  = x(:,16); 
xl_65 = x(:,17); yl_65 = x(:,18);  zl_65 = x(:,19); xl_64 = x(:,20); 
yl_64 = x(:,21); zl_64 = x(:,22); w_5   = x(:,23); xi_5 = x(:,24); 
yi_5  = x(:,25); zi_5  = x(:,26); bi_5  = x(:,27); xl_53 = x(:,28); 
yl_53 = x(:,29); zl_53 = x(:,30); xl_52 = x(:,31); yl_52 = x(:,32); 
zl_52 = x(:,33); w_4   = x(:,53); xi_4  = x(:,54); yi_4  = x(:,55); 
zi_4  = x(:,56); bi_4  = x(:,57); xl_41 = x(:,58); yl_41 = x(:,59); 
zl_41 = x(:,60); xl_42 = x(:,61); yl_42 = x(:,62); zl_42 = x(:,63);
w_3   = x(:,34); xi_3  = x(:,35); yi_3  = x(:,36); zi_3  = x(:,37);
bi_3  = x(:,38); xl_32 = x(:,39); yl_32 = x(:,40); zl_32 = x(:,41); 
xl_3  = x(:,42); yl_3  = x(:,43); zl_3  = x(:,44);
w_2   = x(:,45); xi_2  = x(:,46); yi_2  = x(:,47); zi_2  = x(:,48); 
bi_2  = x(:,49); xl_21 = x(:,50); yl_21 = x(:,51); zl_21 = x(:,52);
w_1 = x(:,64); xi_1 = x(:,65); yi_1 = x(:,66); zi_1 = x(:,67); 
bi_1 = x(:,68); xl_1 = x(:,69); yl_1 = x(:,70); zl_1 = x(:,71);
ca_7 = x(:,72); ca_6 = x(:,76); ca_5 = x(:,80); ca_4 = x(:,84); 
ca_3 = x(:,88); ca_2 = x(:,92); ca_1 = x(:,96);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1)
plot(t, xi_7,'r','LineWidth',3)
hold on
plot(t, yi_7,'b','LineWidth',3)
hold on
plot(t, zi_7,'k', 'LineWidth',3)
hold on
plot(t, xl_75,'--', 'LineWidth',3)
hold on
plot(t, yl_75,'--', 'LineWidth',3)
hold on
plot(t, zl_75,'--','LineWidth',3)
hold on
plot(t, xl_76,'--', 'LineWidth',3)
hold on
plot(t, yl_76,'--', 'LineWidth',3)
hold on
plot(t, zl_76,'--','LineWidth',3)
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('(mM)')
title('Cell_7')
legend('[Na^+]_i','[K^+]_i','[Cl^-]_i', ...
    '[Na^+]_{l_{75}}','[K^+]_{l_{75}}','[Cl^-]_{l_{75}}',...
    '[Na^+]_{l_{76}}','[K^+]_{l_{76}}','[Cl^-]_{l_{76}}')
box off
legend boxoff


figure(2)
plot(t, w_7, 'LineWidth',3)
hold on
plot(t, w_6, 'LineWidth',3)
hold on
plot(t, w_5, 'LineWidth',3)
hold on
plot(t, w_4, 'LineWidth',3)
hold on
plot(t, w_3, 'LineWidth',3)
hold on
plot(t, w_2, 'LineWidth',3)
hold on
plot(t, w_1, 'LineWidth',3)
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('Volume \mu m^3')
title('Cellular Volume')
legend('Cell_7','Cell_6','Cell_5','Cell_4','Cell_3','Cell_2','Cell_1')
box off
legend boxoff

figure(3)
plot(t, ca_7 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_6 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_5 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_4 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_3 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_2 * 1e3, 'LineWidth',3)
hold on
plot(t, ca_1 * 1e3, 'LineWidth',3)
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('[Ca^{2+}]_i (nM)')
legend('Ca^{2+}_7','Ca^{2+}_6','Ca^{2+}_5','Ca^{2+}_4','Ca^{2+}_3','Ca^{2+}_2','Ca^{2+}_1')
box off
legend boxoff

CA = [ca_1(end);ca_2(end);ca_3(end);ca_4(end);ca_5(end);ca_6(end);ca_7(end)]*1e3;

qb_7 = c7.par.Lb * ( 2 * ( xi_7 + yi_7 + c7.par.Hy ) ...
                            + c7.par.CO20 - ( c7.par.Nae + c7.par.Ke ...
                                            + c7.par.Cle + c7.par.HCO3e ) );

qa_76 = c7.par.La * ( 2 * ( xl_76 + yl_76 - xi_7 - yi_7 - c7.par.Hy )...
- c7.par.CO20 + c7.par.Ul );

qt_76 = c7.par.Lt * ( 2 * ( xl_76 + yl_76 ) + c7.par.Ul ...
- ( c7.par.Nae + c7.par.Ke + c7.par.Cle + c7.par.HCO3e ) ); 

qa_75 = c7.par.La * ( 2 * ( xl_75 + yl_75 - xi_7 - yi_7 - c7.par.Hy )...
- c7.par.CO20 + c7.par.Ul );

qt_75 = c7.par.Lt * ( 2 * ( xl_75 + yl_75 ) + c7.par.Ul ...
- ( c7.par.Nae + c7.par.Ke + c7.par.Cle + c7.par.HCO3e ) ); 

qT_76 = qa_76 + qt_76; 

qT_75 = qa_75 + qt_75;

Q_7 = (qT_76 + qT_75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qb_6 = c6.par.Lb * ( 2 * ( xi_6 + yi_6 + c6.par.Hy ) ...
                            + c6.par.CO20 - ( c6.par.Nae + c6.par.Ke ...
                                             + c6.par.Cle + c6.par.HCO3e ) );

qa_67 = c6.par.La * ( 2 * ( xl_76 + yl_76 - xi_6 - yi_6 - c6.par.Hy )...
                                           - c6.par.CO20 + c6.par.Ul );

qt_67 = c6.par.Lt * ( 2 * ( xl_76 + yl_76 ) + c6.par.Ul ...
                    - ( c6.par.Nae + c6.par.Ke + c6.par.Cle + c6.par.HCO3e ) );

qa_65 = c6.par.La * ( 2 * ( xl_65 + yl_65 - xi_6 - yi_6 - c6.par.Hy )...
                                           - c6.par.CO20 + c6.par.Ul );

qt_65 = c6.par.Lt * ( 2 * ( xl_65 + yl_65 ) + c6.par.Ul ...
                    - ( c6.par.Nae + c6.par.Ke + c6.par.Cle + c6.par.HCO3e ) );

qa_64 = c6.par.La * ( 2 * ( xl_64 + yl_64 - xi_6 - yi_6 - c6.par.Hy )...
                                           - c6.par.CO20 + c6.par.Ul );

qt_64 = c6.par.Lt * ( 2 * ( xl_64 + yl_64 ) + c6.par.Ul ...
                    - ( c6.par.Nae + c6.par.Ke + c6.par.Cle + c6.par.HCO3e ) );


qT_67 = qa_67 + qt_67; 

qT_65 = qa_65 + qt_65; 

qT_64 = qa_64 + qt_64;

Q_6 = (qT_67 + qT_65 + qT_64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qb_5 = c5.par.Lb * ( 2 * ( xi_5 + yi_5 + c5.par.Hy ) ...
                               + c5.par.CO20 - ( c5.par.Nae + c5.par.Ke ...
                                             + c5.par.Cle + c5.par.HCO3e ) );

qa_57 = c5.par.La * ( 2 * ( xl_75 + yl_75 - xi_5 - yi_5 - c5.par.Hy )...
                                           - c5.par.CO20 + c5.par.Ul );

qt_57 = c5.par.Lt * ( 2 * ( xl_75 + yl_75 ) + c5.par.Ul ...
                    - ( c5.par.Nae + c5.par.Ke + c5.par.Cle + c5.par.HCO3e ) );

qa_56 = c5.par.La * ( 2 * ( xl_65 + yl_65 - xi_5 - yi_5 - c5.par.Hy )...
                                           - c5.par.CO20 + c5.par.Ul );

qt_56 = c5.par.Lt * ( 2 * ( xl_65 + yl_65 ) + c5.par.Ul ...
                    - ( c5.par.Nae + c5.par.Ke + c5.par.Cle + c5.par.HCO3e ) );

qa_53 = c5.par.La * ( 2 * ( xl_53 + yl_53 - xi_5 - yi_5 - c5.par.Hy )...
                                           - c5.par.CO20 + c5.par.Ul );

qt_53 = c5.par.Lt * ( 2 * ( xl_53 + yl_53 ) + c5.par.Ul ...
                    - ( c5.par.Nae + c5.par.Ke + c5.par.Cle + c5.par.HCO3e ) );

qa_52 = c5.par.La * ( 2 * ( xl_52 + yl_52 - xi_5 - yi_5 - c5.par.Hy )...
                                           - c5.par.CO20 + c5.par.Ul );

qt_52 = c5.par.Lt * ( 2 * ( xl_52 + yl_52 ) + c5.par.Ul ...
                    - ( c5.par.Nae + c5.par.Ke + c5.par.Cle + c5.par.HCO3e ) );

qT_57 = qa_57 + qt_57; 
qT_56 = qa_56 + qt_56; 
qT_53 = qa_53 + qt_53;
qT_52 = qa_52 + qt_52;
Q_5 = ( qT_57 + qT_56 + qT_53 + qT_52 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qb_4 = c4.par.Lb * ( 2 * ( xi_4 + yi_4 + c4.par.Hy ) ...
                            + c4.par.CO20 - ( c4.par.Nae + c4.par.Ke ...
                                             + c4.par.Cle + c4.par.HCO3e ) );

qa_46 = c4.par.La * ( 2 * ( xl_64 + yl_64 - xi_4 - yi_4 - c4.par.Hy )...
                                           - c4.par.CO20 + c4.par.Ul );

qt_46 = c4.par.Lt * ( 2 * ( xl_64 + yl_64 ) + c4.par.Ul ...
                    - ( c4.par.Nae + c4.par.Ke + c4.par.Cle + c4.par.HCO3e ) );

qa_42 = c4.par.La * ( 2 * ( xl_42 + yl_42 - xi_4 - yi_4 - c4.par.Hy )...
                                           - c4.par.CO20 + c4.par.Ul );

qt_42 = c4.par.Lt * ( 2 * ( xl_42 + yl_42 ) + c4.par.Ul ...
                    - ( c4.par.Nae + c4.par.Ke + c4.par.Cle + c4.par.HCO3e ) );

qa_41 = c4.par.La * ( 2 * ( xl_41 + yl_41 - xi_4 - yi_4 - c4.par.Hy )...
                                           - c4.par.CO20 + c4.par.Ul );

qt_41 = c4.par.Lt * ( 2 * ( xl_41 + yl_41 ) + c4.par.Ul ...
                    - ( c4.par.Nae + c4.par.Ke + c4.par.Cle + c4.par.HCO3e ) );


qT_46 = qa_46 + qt_46; 

qT_42 = qa_42 + qt_42; 

qT_41 = qa_41 + qt_41;

Q_4 = ( qT_46 + qT_42 + qT_41 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


qb_3 = c3.par.Lb * ( 2 * ( xi_3 + yi_3 + c3.par.Hy ) ...
+ c3.par.CO20 - ( c3.par.Nae + c3.par.Ke ...
+ c3.par.Cle + c3.par.HCO3e ) );

qa_3 = c3.par.La * ( 2 * ( xl_3 + yl_3 - xi_3 - yi_3 - c3.par.Hy )...
- c3.par.CO20 + c3.par.Ul );

qt_3 = c3.par.Lt * ( 2 * ( xl_3 + yl_3 ) + c3.par.Ul ...
- ( c3.par.Nae + c3.par.Ke + c3.par.Cle + c3.par.HCO3e ) );

qa_35 = c3.par.La * ( 2 * ( xl_53 + yl_53 - xi_3 - yi_3 - c3.par.Hy )...
- c3.par.CO20 + c3.par.Ul );

qt_35 = c3.par.Lt * ( 2 * ( xl_53 + yl_53 ) + c3.par.Ul ...
- ( c3.par.Nae + c3.par.Ke + c3.par.Cle + c3.par.HCO3e ) );


qa_32 = c3.par.La * ( 2 * ( xl_32 + yl_32 - xi_3 - yi_3 - c3.par.Hy )...
- c3.par.CO20 + c3.par.Ul );

qt_32 = c3.par.Lt * ( 2 * ( xl_32 + yl_32 ) + c3.par.Ul ...
- ( c3.par.Nae + c3.par.Ke + c3.par.Cle + c3.par.HCO3e ) );

qT_3 = qa_3 + qt_3; 

qT_35 = qa_35 + qt_35; 

qT_32 = qa_32 + qt_32;

Q_3 = ( qT_35 + qT_32 + qT_3 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qb_2 = c2.par.Lb * ( 2 * ( xi_2 + yi_2 + c2.par.Hy ) ...
+ c2.par.CO20 - ( c2.par.Nae + c2.par.Ke ...
+ c2.par.Cle + c2.par.HCO3e ) );

qa_25 = c2.par.La * ( 2 * ( xl_52 + yl_52 - xi_2 - yi_2 - c2.par.Hy )...
- c2.par.CO20 + c2.par.Ul );

qt_25 = c2.par.Lt * ( 2 * ( xl_52 + yl_52 ) + c2.par.Ul ...
- ( c2.par.Nae + c2.par.Ke + c2.par.Cle + c2.par.HCO3e ) );

qa_24 = c2.par.La * ( 2 * ( xl_42 + yl_42 - xi_2 - yi_2 - c2.par.Hy )...
- c2.par.CO20 + c2.par.Ul );

qt_24 = c2.par.Lt * ( 2 * ( xl_42 + yl_42 ) + c2.par.Ul ...
- ( c2.par.Nae + c2.par.Ke + c2.par.Cle + c2.par.HCO3e ) );

qa_23 = c2.par.La * ( 2 * ( xl_32 + yl_32 - xi_2 - yi_2 - c2.par.Hy )...
- c2.par.CO20 + c2.par.Ul );

qt_23 = c2.par.Lt * ( 2 * ( xl_32 + yl_32 ) + c2.par.Ul ...
- ( c2.par.Nae + c2.par.Ke + c2.par.Cle + c2.par.HCO3e ) );

qa_21 = c2.par.La * ( 2 * ( xl_21 + yl_21 - xi_2 - yi_2 - c2.par.Hy )...
- c2.par.CO20 + c2.par.Ul );

qt_21 = c2.par.Lt * ( 2 * ( xl_21 + yl_21 ) + c2.par.Ul ...
- ( c2.par.Nae + c2.par.Ke + c2.par.Cle + c2.par.HCO3e ) );


qT_25 = qa_25 + qt_25; 

qT_24 = qa_24 + qt_24; 

qT_23 = qa_23 + qt_23;

qT_21 = qa_21 + qt_21; 

Q_2 = ( qT_25 + qT_24 + qT_23 + qT_21 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qb_1 = c1.par.Lb * ( 2 * ( xi_1 + yi_1 + c1.par.Hy ) ...
+ c1.par.CO20 - ( c1.par.Nae + c1.par.Ke ...
+ c1.par.Cle + c1.par.HCO3e ) );

qa_1 = c1.par.La * ( 2 * ( xl_1 + yl_1 - xi_1 - yi_1 - c1.par.Hy )...
- c1.par.CO20 + c1.par.Ul );

qt_1 = c1.par.Lt * ( 2 * ( xl_1 + yl_1 ) + c1.par.Ul ...
- ( c1.par.Nae + c1.par.Ke + c1.par.Cle + c1.par.HCO3e ) ); 

qa_14 = c1.par.La * ( 2 * ( xl_41 + yl_41 - xi_1 - yi_1 - c1.par.Hy )...
- c1.par.CO20 + c1.par.Ul );

qt_14 = c1.par.Lt * ( 2 * ( xl_41 + yl_41 ) + c1.par.Ul ...
- ( c1.par.Nae + c1.par.Ke + c1.par.Cle + c1.par.HCO3e ) );

qa_12 = c1.par.La * ( 2 * ( xl_21 + yl_21 - xi_1 - yi_1 - c1.par.Hy )...
- c1.par.CO20 + c1.par.Ul );

qt_12 = c1.par.Lt * ( 2 * ( xl_21 + yl_21 ) + c1.par.Ul ...
- ( c1.par.Nae + c1.par.Ke + c1.par.Cle + c1.par.HCO3e ) );

qT_1 = ( qa_1 + qt_1 );

qT_14 = ( qa_14 + qt_14 );
 
qT_12 = ( qa_12 + qt_12 );

Q_1 = ( qT_1 + qT_12 + qT_14 );

% figure()
% plot(t, Q_7+3, 'LineWidth',3)
% hold on
% plot(t, Q_6+3, 'LineWidth',3)
% hold on
% plot(t, Q_5+3, 'LineWidth',3)
% hold on
% plot(t, Q_4+3, 'LineWidth',3)
% hold on
% plot(t, Q_3+3, 'LineWidth',3)
% hold on
% plot(t, Q_2+3, 'LineWidth',3)
% hold on
% plot(t, Q_1+3, 'LineWidth',3)
% ax = gca;
% ax.LineWidth = 2.1;
% ax.FontSize = 20;
% xlabel('Time (sec)')
% ylabel('Individual FFR (\mu m^3 /sec)')
% legend('Cell_7','Cell_6','Cell_5','Cell_4','Cell_3','Cell_2','Cell_1')
% box off
% legend boxoff
% saveas(gcf,'Cell_FF','epsc')

Q_76 = (qT_76 + qT_67);
Q_75 = (qT_75 + qT_57);
Q_65 = (Q_76 + Q_75 + qT_65 + qT_56); 
Q_53 = (qT_53 + qT_35);
Q_52 = (Q_65 + qT_52 + qT_25);
Q_21 = (qT_21 + qT_12 + qT_1);
Q_41 = (qT_41 + qT_14);
Q_64 = (qT_64 + qT_46);
Q_42 = qT_42 + qT_24 + Q_64 + Q_41 + Q_21;
Q_32 = (Q_52 + Q_53 + qT_32 + qT_23 + Q_42);
QT_cluster = (Q_32 + qT_3) - 20;


Time_cluster = t;
figure()
plot(t, QT_cluster, 'LineWidth',3)
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('Cluster FFR (\mu m^3 /sec)')
box off