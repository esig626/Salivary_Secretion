classdef Cluster
   
    properties 
    Q
    j
    q
    intra
    luminal
    Eq
    end
    % Constructor Method for Cluster
    methods
        function cluster = Cluster(c)
        [cluster.Q,cluster.j,cluster.q,cluster.luminal] ...
                                  = qj(c(1),c(2),c(3),c(4),c(5),c(6),c(7));  
        for i = 1:7
            cluster.intra{i,1} = Cluster.intracell(c(i));
        end
        
        cluster.Eq = equations(cluster.intra,cluster.luminal);
        
        function [q,J,qq,lum] = qj(c1,c2,c3,c4,c5,c6,c7)
            q.f76 = (c7.fx.qT(1) + c6.fx.qT(1));
            q.f75 = (c7.fx.qT(2) + c5.fx.qT(1));
            q.f65 = (q.f76 + q.f75 + c6.fx.qT(2) + c5.fx.qT(2));
            q.f53 = (c5.fx.qT(3) + c3.fx.qT(1));
            q.f52 = (q.f65 + c5.fx.qT(4) + c2.fx.qT(1));
            q.f21 = (c2.fx.qT(4) + c1.fx.qT(2) + c1.fx.qT(3));
            q.f41 = (c4.fx.qT(3) + c1.fx.qT(1));
            q.f64 = (c6.fx.qT(3) + c4.fx.qT(1));
            q.f42 = c4.fx.qT(2) + c2.fx.qT(2) + q.f64 + q.f41 + q.f21;
            q.f32 = (q.f52 + q.f53 + c3.fx.qT(3) + c2.fx.qT(3) + q.f42);
            q.f33 = (q.f32 + c3.fx.qT(2));
            
            J.jtx =[ c1.fx.jtx(3), 0;
                (c2.fx.jtx(4) + c1.fx.jtx(2)), 0;
                (c3.fx.jtx(3) + c2.fx.jtx(3)),c3.fx.jtx(2);
                (c4.fx.jtx(3) + c1.fx.jtx(1)),(c4.fx.jtx(2) + c2.fx.jtx(2));
                (c5.fx.jtx(3) + c3.fx.jtx(1)),(c5.fx.jtx(4) + c2.fx.jtx(1));
                (c6.fx.jtx(2) + c5.fx.jtx(2)),(c6.fx.jtx(3) + c4.fx.jtx(1));
                (c7.fx.jtx(1) + c6.fx.jtx(1)),(c7.fx.jtx(2) + c5.fx.jtx(1))];
            
            J.jty = [ c1.fx.jty(3), 0;
                (c2.fx.jty(4) + c1.fx.jty(2)), 0;
                (c3.fx.jty(3) + c2.fx.jty(3)),c3.fx.jty(2);
                (c4.fx.jty(3) + c1.fx.jty(1)),(c4.fx.jty(2) + c2.fx.jty(2));
                (c5.fx.jty(3) + c3.fx.jty(1)),(c5.fx.jty(4) + c2.fx.jty(1));
                (c6.fx.jty(2) + c5.fx.jty(2)),(c6.fx.jty(3) + c4.fx.jty(1));
                (c7.fx.jty(1) + c6.fx.jty(1)),(c7.fx.jty(2) + c5.fx.jty(1))];
            
            J.jz =[- c1.fx.jz(3), 0;
                -(c2.fx.jz(4) + c1.fx.jz(2)), 0;
                -(c3.fx.jz(3) + c2.fx.jz(3)),-c3.fx.jz(2);
                -(c4.fx.jz(3) + c1.fx.jz(1)),-(c4.fx.jz(2)  + c2.fx.jz(2));
                -(c5.fx.jz(3) + c3.fx.jz(1)),-(c5.fx.jz(4) +  c2.fx.jz(1));
                -(c6.fx.jz(2) + c5.fx.jz(2)),-(c6.fx.jz(3) + c4.fx.jz(1));
                -(c7.fx.jz(1) + c6.fx.jz(1)),-(c7.fx.jz(2) + c5.fx.jz(1))];
            
            qq.qx =[-c1.fx.qT(3) * c1.var.xl(3), 0;
                c1.fx.qT(3) * c1.var.xl(3) - q.f21 * c2.var.xl(3), 0;
                q.f52 * c5.var.xl(4) + q.f53 * c5.var.xl(3) + q.f42 * c4.var.xl(2) - q.f32 * c3.var.xl(3), q.f32 * c3.var.xl(3) - q.f33 * c3.var.xl(2);
                -q.f41 * c4.var.xl(3),+ q.f64 * c6.var.xl(3) + q.f41 * c4.var.xl(3) + q.f21 * c2.var.xl(3) - q.f42 * c4.var.xl(2);
                -q.f53 * c5.var.xl(3), q.f65 * c6.var.xl(2) - q.f52 * c5.var.xl(4);
                q.f76 * c7.var.xl(1) + q.f75 * c7.var.xl(2) - q.f65 * c6.var.xl(2),- q.f64 * c6.var.xl(3);
                -q.f76 * c7.var.xl(1),- q.f75 * c7.var.xl(2)];
            
            qq.qy =[-c1.fx.qT(3) * c1.var.yl(3), 0;
                c1.fx.qT(3) * c1.var.yl(3) - q.f21 * c2.var.yl(3), 0;
                q.f52 * c5.var.yl(4) + q.f53 * c5.var.yl(3) + q.f42 * c4.var.yl(2) - q.f32 * c3.var.yl(3), q.f32 * c3.var.yl(3) - q.f33 * c3.var.yl(2);
                -q.f41 * c4.var.yl(3),q.f64 * c6.var.yl(3) + q.f41 * c4.var.yl(3) + q.f21 * c2.var.yl(3) - q.f42 * c4.var.yl(2);
                -q.f53 * c5.var.yl(3),q.f65 * c6.var.yl(2) - q.f52 * c5.var.yl(4);
                q.f76 * c7.var.yl(1) + q.f75 * c7.var.yl(2) - q.f65 * c6.var.yl(2),- q.f64 * c6.var.yl(3);
                -q.f76 * c7.var.yl(1),-q.f75 * c7.var.yl(2)];
            
            qq.qz =[-c1.fx.qT(3) * c1.var.zl(3), 0;
                c1.fx.qT(3) * c1.var.zl(3) - q.f21 * c2.var.zl(3),0;
                q.f52 * c5.var.zl(4) + q.f53 * c5.var.zl(3) + q.f42 * c4.var.zl(2) - q.f32 * c3.var.zl(3), q.f32 * c3.var.zl(3) - q.f33 * c3.var.zl(2);
                -q.f41 * c4.var.zl(3), q.f64 * c6.var.zl(3) + q.f41 * c4.var.zl(3) + q.f21 * c2.var.zl(3) - q.f42 * c4.var.zl(2);
                -q.f53 * c5.var.zl(3),+ q.f65 * c6.var.zl(2) - q.f52 * c5.var.zl(4);
                q.f76 * c7.var.zl(1) + q.f75 * c7.var.zl(2) - q.f65 * c6.var.zl(2),- q.f64 * c6.var.zl(3);
                -q.f76 * c7.var.zl(1), -q.f75 * c7.var.zl(2)];
            
            lum.dxl = J.jtx + qq.qx;
            lum.dyl = J.jty + qq.qy;
            lum.dzl = J.jz  + qq.qz;
        end
            function dx = equations(intra,luminal)
                dx(6) = luminal.dxl(7,1);
                dx(9) = luminal.dxl(7,2);
                dx(7)  = luminal.dyl(7,1);
                dx(10) = luminal.dyl(7,2);
                dx(8)  = luminal.dzl(7,1);
                dx(11) = luminal.dzl(7,2);
                dx(17) = luminal.dxl(6,1);
                dx(20) = luminal.dxl(6,2);
                dx(18) = luminal.dyl(6,1);
                dx(21) = luminal.dyl(6,2);
                dx(19) = luminal.dzl(6,1);
                dx(22) = luminal.dzl(6,2);
                dx(28) = luminal.dxl(5,1);
                dx(31) = luminal.dxl(5,2);
                dx(29) = luminal.dyl(5,1);
                dx(32) = luminal.dyl(5,2);
                dx(30) = luminal.dzl(5,1);
                dx(33) = luminal.dzl(5,2);
                dx(58) = luminal.dxl(4,1);
                dx(61) = luminal.dxl(4,2);
                dx(59) = luminal.dyl(4,1);
                dx(62) = luminal.dyl(4,2);
                dx(60) = luminal.dzl(4,1);
                dx(63) = luminal.dzl(4,2);
                dx(39) = luminal.dxl(3,1);
                dx(42) = luminal.dxl(3,2);
                dx(40) = luminal.dyl(3,1);
                dx(43) = luminal.dyl(3,2);
                dx(41) = luminal.dzl(3,1);
                dx(44) = luminal.dzl(3,2);
                dx(50) = luminal.dxl(2,1);
                dx(51) = luminal.dyl(2,1);
                dx(52) = luminal.dzl(2,1);
                dx(69) = luminal.dxl(1,1);
                dx(70) = luminal.dyl(1,1);
                dx(71) = luminal.dzl(1,1);
                dx(1) = intra{7}.dw;
                dx(2) = intra{7}.dx;
                dx(3) = intra{7}.dy;
                dx(4) = intra{7}.dz;
                dx(5) = intra{7}.db;
                dx(12) = intra{6}.dw;
                dx(13) = intra{6}.dx;
                dx(14) = intra{6}.dy;
                dx(15) = intra{6}.dz;
                dx(16) = intra{6}.db;
                dx(23) = intra{5}.dw;
                dx(24) = intra{5}.dx;
                dx(25) = intra{5}.dy;
                dx(26) = intra{5}.dz;
                dx(27) = intra{5}.db;
                dx(53) = intra{4}.dw;
                dx(54) = intra{4}.dx;
                dx(55) = intra{4}.dy;
                dx(56) = intra{4}.dz;
                dx(57) = intra{4}.db;
                dx(34) = intra{3}.dw;
                dx(35) = intra{3}.dx;
                dx(36) = intra{3}.dy;
                dx(37) = intra{3}.dz;
                dx(38) = intra{3}.db;
                dx(45) = intra{2}.dw;
                dx(46) = intra{2}.dx;
                dx(47) = intra{2}.dy;
                dx(48) = intra{2}.dz;
                dx(49) = intra{2}.db;
                dx(64) = intra{1}.dw;
                dx(65) = intra{1}.dx;
                dx(66) = intra{1}.dy;
                dx(67) = intra{1}.dz;
                dx(68) = intra{1}.db;
                dx(72) = intra{7}.dca;
                dx(73) = intra{7}.dip;
                dx(74) = intra{7}.dhh;
                dx(75) = intra{7}.dgg;
                dx(76) = intra{6}.dca;
                dx(77) = intra{6}.dip;
                dx(78) = intra{6}.dhh;
                dx(79) = intra{6}.dgg;
                dx(80) = intra{5}.dca;
                dx(81) = intra{5}.dip;
                dx(82) = intra{5}.dhh;
                dx(83) = intra{5}.dgg;
                dx(84) = intra{4}.dca;
                dx(85) = intra{4}.dip;
                dx(86) = intra{4}.dhh;
                dx(87) = intra{4}.dgg;
                dx(88) = intra{3}.dca;
                dx(89) = intra{3}.dip;
                dx(90) = intra{3}.dhh;
                dx(91) = intra{3}.dgg;
                dx(92) = intra{2}.dca;
                dx(93) = intra{2}.dip;
                dx(94) = intra{2}.dhh;
                dx(95) = intra{2}.dgg;
                dx(96) = intra{1}.dca;
                dx(97) = intra{1}.dip;
                dx(98) = intra{1}.dhh;
                dx(99) = intra{1}.dgg;
            end
        end
        
    end
    methods (Static)
        function in = intracell(c)
            in.dw = c.fx.jw;
            in.dx =(c.fx.jNk1-3* c.fx.jNaK+c.fx.j1-c.fx.j4-c.fx.jw*c.var.xi)/c.var.w ;
            in.dy =(c.fx.jNk1+2*c.fx.jNaK-c.fx.jy-c.fx.jw * c.var.yi)/c.var.w;
            in.dz =(2 * c.fx.jNk1+c.fx.j4+sum(c.fx.jz)-c.fx.jw*c.var.zi)/c.var.w;
            in.db =(c.fx.jB-2*c.fx.j4-c.fx.jw*c.var.bi)/c.var.w;
            in.dca = c.fx.jcalcium;
            in.dip = c.fx.jip3;
            in.dhh = c.fx.jhh;
            in.dgg = c.fx.jgg;          
        end
    end
        
end






































































