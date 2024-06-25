function dx = Cluster_Secretion(t,x,c7,c6,c5,c4,c3,c2,c1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t
% CALCIUM MODEL ODE (Closed)
ca_7 = x(72); ip_7 = x(73); hh_7 = x(74); gg_7 = x(75); 
ca_6 = x(76); ip_6 = x(77); hh_6 = x(78); gg_6 = x(79);
ca_5 = x(80); ip_5 = x(81); hh_5 = x(82); gg_5 = x(83);
ca_4 = x(84); ip_4 = x(85); hh_4 = x(86); gg_4 = x(87);
ca_3 = x(88); ip_3 = x(89); hh_3 = x(90); gg_3 = x(91);
ca_2 = x(92); ip_2 = x(93); hh_2 = x(94); gg_2 = x(95);
ca_1 = x(96); ip_1 = x(97); hh_1 = x(98); gg_1 = x(99);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cell 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_7   = x(1); xi_7  = x(2); yi_7  = x(3); zi_7  = x(4); bi_7  = x(5); xl_76 = x(6); yl_76 = x(7); zl_76 = x(8); xl_75 = x(9); yl_75 = x(10); zl_75 = x(11);
%%%%% Cell 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_6   = x(12); xi_6  = x(13); yi_6  = x(14); zi_6  = x(15); bi_6  = x(16); xl_65 = x(17); yl_65 = x(18); zl_65 = x(19); xl_64 = x(20); yl_64 = x(21); zl_64 = x(22); 
%%%%% Cell 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_5   = x(23); xi_5  = x(24); yi_5  = x(25); zi_5  = x(26); bi_5  = x(27); xl_53 = x(28); yl_53 = x(29); zl_53 = x(30); xl_52 = x(31); yl_52 = x(32); zl_52 = x(33);
%%%%% Cell 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_4   = x(53); xi_4  = x(54); yi_4  = x(55); zi_4  = x(56); bi_4  = x(57); xl_41 = x(58); yl_41 = x(59); zl_41 = x(60); xl_42 = x(61); yl_42 = x(62); zl_42 = x(63); 
%%%%% Cell 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_3   = x(34); xi_3  = x(35); yi_3  = x(36); zi_3  = x(37); bi_3  = x(38); xl_32 = x(39); yl_32 = x(40); zl_32 = x(41); xl_3  = x(42); yl_3  = x(43); zl_3  = x(44);
%%%% Cell 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_2   = x(45); xi_2  = x(46); yi_2  = x(47); zi_2  = x(48); bi_2  = x(49); xl_21 = x(50); yl_21 = x(51); zl_21 = x(52);
%%%% Cell 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_1 = x(64); xi_1 = x(65); yi_1 = x(66); zi_1 = x(67); bi_1 = x(68); xl_1 = x(69); yl_1 = x(70); zl_1 = x(71);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cell_7 = [w_7, xi_7, yi_7, zi_7, bi_7, ca_7, ip_7, hh_7, gg_7, xl_76, yl_76, zl_76, xl_75, yl_75, zl_75];
Cell_6 = [w_6, xi_6, yi_6, zi_6, bi_6, ca_6, ip_6, hh_6, gg_6, xl_76, yl_76, zl_76, xl_65, yl_65, zl_65, xl_64, yl_64, zl_64];
Cell_5 = [w_5, xi_5, yi_5, zi_5, bi_5, ca_5, ip_5, hh_5, gg_5, xl_75, yl_75, zl_75, xl_65, yl_65, zl_65, xl_53, yl_53, zl_53, xl_52, yl_52, zl_52];
Cell_4 = [w_4, xi_4, yi_4, zi_4, bi_4, ca_4, ip_4, hh_4, gg_4, xl_64, yl_64, zl_64, xl_42, yl_42, zl_42, xl_41, yl_41, zl_41];
Cell_3 = [w_3, xi_3, yi_3, zi_3, bi_3, ca_3, ip_3, hh_3, gg_3, xl_53, yl_53, zl_53, xl_3, yl_3, zl_3, xl_32, yl_32, zl_32];
Cell_2 = [w_2, xi_2, yi_2, zi_2, bi_2, ca_2, ip_2, hh_2, gg_2, xl_52, yl_52, zl_52, xl_42, yl_42, zl_42, xl_32, yl_32, zl_32, xl_21, yl_21, zl_21];
Cell_1 = [w_1, xi_1, yi_1, zi_1, bi_1, ca_1, ip_1, hh_1, gg_1, xl_41, yl_41, zl_41, xl_21, yl_21, zl_21, xl_1, yl_1, zl_1];


[c7.var,c7.fx] = Cell.fluxes(c7.c_no,c7.ngh,Cell_7,c7.par,1,1);

[c6.var,c6.fx] = Cell.fluxes(c6.c_no,c6.ngh,Cell_6,c6.par,...
                                                c7.fx.vtx(1),c7.fx.vty(1));
[c5.var,c5.fx] = Cell.fluxes(c5.c_no,c5.ngh,Cell_5,c5.par,...
                [c7.fx.vtx(2);c6.fx.vtx(2)],[c7.fx.vty(2);c6.fx.vty(2)]);
[c4.var,c4.fx] = Cell.fluxes(c4.c_no,c4.ngh,Cell_4,c4.par,...
                                                c6.fx.vtx(3),c6.fx.vty(3));
[c3.var,c3.fx] = Cell.fluxes(c3.c_no,c3.ngh,Cell_3,c3.par,...
                                                c5.fx.vtx(3),c5.fx.vty(3));
[c2.var,c2.fx] = Cell.fluxes(c2.c_no,c2.ngh,Cell_2,c2.par,...
                        [c5.fx.vtx(4);c4.fx.vtx(2);c3.fx.vtx(3)],...
                                 [c5.fx.vty(4);c4.fx.vty(2);c3.fx.vty(3)]);
[c1.var,c1.fx] = Cell.fluxes(c1.c_no,c1.ngh,Cell_1,c1.par,...
                                        [c4.fx.vtx(3);c2.fx.vtx(4)],...
                                            [c4.fx.vty(3);c2.fx.vty(4)]);
                                                                            
c=[c1,c2,c3,c4,c5,c6,c7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl = Cluster(c);

dx=cl.Eq';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end