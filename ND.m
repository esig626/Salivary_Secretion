function F = ND(t,U,V_PLC)
% if t<60
%     V_PLC = 0;
% end
% ---------------------------------------------------
% DECLARE VARIABLES
% ---------------------------------------------------
    V   = U(1);            
    XI  = U(2);           
    YI  = U(3);            
    ZI  = U(4);           
    BI  = U(5);                     
    XL  = U(6);           
    YL  = U(7);
    CA  = U(8);
    IP  = U(9);
    H   = U(10);
    G   = U(11);
    CT  = 3/V;U(12);
%     HI  = U(12);

%     0.0715911332507813;
% ---------------------------------------------------
% CALCIUM PARAMETERS
% ---------------------------------------------------
    K_F = 45;
    K_BETA =0.4;
    K_P = 0.2;
    K_C = 0.2;
    K_H = 0.08;
    V_P = 0.9;
    K_BAR = 0.00001957;
    GAMMA=5.5;
    K_PLC = 0.07;
    K_3K=0.4;
    V_3K = 0.05;
    V_5K = 0.05;
    TAU_MAX = 100;
    K_TAU = 0.1;
    V_RYR = 0.07;
    K_RYR = 0.2;
    K_RYR2 = 0.15;
    TAU_RYR = 1;
    ALPHA_0 = 0.0027;
    ALPHA_1 = 0.07;
    ALPHA_2 = 15;
    DELTA = 0;
    V_PM = 0.11;
    K_PM = 0.3;
    K_CE = 8;
    SCALE_HY = 10^2;
% ---------------------------------------------------
% FF MODEL PARAMETERS
% ---------------------------------------------------
    ANAK = 0.1606;
    ANKCC1 = 2.0712;
    GTNA = 2.0179E+06;
    GTK = 1.4600E+05;
    GCL = 1.8239E+09;
    GK = 1.6722E+08;
    G1 = 90.8356;
    G4 = 1.3361E+03;
    GB = 9.7680;
    ST = 57.6215;
    SB = 324.5745;
    SA = 54.3688;
    V0 = 904;
    KCACC =0.26;     
    KCAKC =0.3;      
    ETA1  =4.49;      
    ETA2  =2;     
    R = 0.001305;       
    ALPHA1 = 0.000641;  
    FAR = 96485.3365;    
    RTF = 26.7137302360965;
    KN = 316.8;
    KP = 0.132;
    A1 = 157.55;                 
    A2 = 2.0096E-05;   
    A3 = 1.0306;                      
    A4 = 1.3852E-06;     
    CLE   =102.6;         
    NAE   =140.2;        
    KE    =5.3;             
    HEL   =3.8904514499428E-05;     
    HCO3E =40;
    KCL = 5.6;
    KB = 100;
    KNA = 15;
    KH = 4.5E-4;
    LA = 4.43354480372293;
    LB = 1.94580852061638;
    LT = 0.0593663800587253;
    UL = 51.7310374376987;   
    CO20 = 6.60012219974207;   
    QA0  = 17.1090179363441;
    VL     = 0.02;
% ---------------------------------------------------
% FUNCTIONS AND NON DIMENSIONAL PARAMETER DEFINITIONS
% ---------------------------------------------------
    NU     = 1/V0; 
    ZL     = XL + YL;
    E_I    = SB * G1 / QA0;
    E_L    = SA * 20 * G1 / QA0;
    XE     = NAE/E_I;
    YE     = KE/E_I;
    ZE     = CLE/E_I;
    HE     = HEL/E_I;
    BE     = HCO3E/E_I;
    C2I    = CO20/E_I;
    PL     = UL/E_L;
    CE     = (CT-CA)*GAMMA;
% ---------------------------------------------------
% CALCIUM MODEL FUNCTIONS
% ---------------------------------------------------
    PHI_C=CA^4/(CA^4+K_C^4);
    PHI_P=IP^2/(K_P^2+IP^2);
    PHI_P_DOWN=K_P^2/(K_P^2+IP^2);
    H_INF=K_H^4/(K_H^4+CA^4);
    TAU = TAU_MAX*(K_TAU^4)/(K_TAU^4+CA^4);
    BETA = PHI_P*PHI_C*H;
    ALPHA = PHI_P_DOWN*(1-PHI_C*H_INF);
    PO=BETA/(BETA + K_BETA*(BETA+ALPHA));
    J_IPR = K_F*PO;
    J_SERCA=V_P*(CA^2-K_BAR*(CE^2))/(CA^2+K_P^2);
    J_RYR = V_RYR *(CA.^4/(CA.^4 + K_RYR.^4))*G;
    G_INF = K_RYR2^2/(K_RYR2^2 + CA^2);
    J_IN = ALPHA_0 + ALPHA_1*(K_CE)^4./((K_CE)^4+CE.^4) + ALPHA_2*V_PLC;
    J_PM = V_PM*CA.^2./((K_PM)^2+CA.^2);
    IP_DEG = (V_5K+V_3K*((CA^2)/(CA^2+K_3K^2)))*IP;
    IP_PLC = V_PLC*((CA^2)/(CA^2+K_PLC^2));
    
% ---------------------------------------------------
% FF MODEL FUNCTIONS
% ---------------------------------------------------
    XV = 46.3677365879351;46.4496200426372;46.2653795533866;
% 
    HI = 0.00000715911332507813;- (XI + YI - ZI - BI - XV./(V*904)).*SCALE_HY;
    LAMBDA_I = 1/(V*V0*E_I);
    LAMBDA_L = 1/(VL*V0*E_L);
    %LAMBDA_IH = 10^6 / (V * V0 * E_I);
    PZ    = 1/(1+(KCACC/CA)^ETA1);
    PY    = 1/(1+(KCAKC/CA)^ETA2);
    KNAKA = R * SB * ANAK * (E_I*E_I*E_I);
    KNAKB = ALPHA1 * E_I;
    JNAK  = KNAKA * (XI*XI*XI) * YE^2 / (YE^2 + KNAKB * (XI*XI*XI));
    JNK1  = SB * ANKCC1 * (A1 - A2 * XI * YI * ZI^2 * E_I^4)...
         /(A3 + A4 * XI * YI * ZI^2 * E_I^4);
    K4  = SB * G4 * E_I^4;
    J4  = (K4 * BI^2 * XI * ZE)/((KB + BI*E_I)^2*(KNA + E_I*XI)*(KCL + E_I*ZE));
    KPI = V0 * GB * KP * E_I;
    KNB = V0 * GB * KN * E_I^2;
    JB  = V *(C2I * KPI - BI * HI /SCALE_HY * KNB);
    N1P = SB * G1 * E_I^2 * XE;
    N1M = SB * G1 * E_I^2 * HE;
    J1  = ((HI /SCALE_HY * N1P)/((KH + HI /SCALE_HY * E_I)*(KNA + XE * E_I))...
       - (N1M * XI)/((KH + HE * E_I)*(KNA + XI * E_I)));
    QB  = LB*(2*((XI * E_I)+(YI * E_I)+(HI /SCALE_HY * E_I))...
       +(C2I * E_I) -((XE * E_I)+(YE * E_I)+(ZE * E_I)+(BE * E_I))) ;            
    QA  = LA*(2*((XL * E_L)+(YL * E_L)-(XI * E_I)...
       -(YI * E_I) -(HI /SCALE_HY * E_I))-(C2I * E_I)+(PL * E_L));
    QTT = LT*(2*((XL * E_L)+(YL * E_L))+(PL * E_L)...
       -((XE * E_I)+(YE * E_I)+(ZE * E_I)+(BE * E_I)));
    JW  = QB - QA;
    QT  = QA + QTT;
    VZ  = RTF * log((ZL * E_L ) / ( ZI * E_I ) );
    VY  = RTF * log((YE * E_I ) / ( YI * E_I ) );
    VTX = RTF * log((XL * E_L) / ( XE * E_I ) );
    VTY = RTF * log((YL * E_L) / ( YE * E_I ) );
    VA  = -(JNAK*FAR*ST*GTK + JNAK*FAR*ST*GTNA ...
       + PY*PZ*GK*GCL*VZ - PY*GK*ST*GTNA*VTX ...
       - PY*GK*ST*GTK*VY - PY*GK*ST*GTNA*VY ...
       - PY*GK*ST*GTK*VTY + PZ*ST*GCL*GTK*VZ ...
       + PZ*ST*GCL*GTNA*VZ) /(PY*GK*ST*GTK ...
       + PY*GK*ST*GTNA + PZ*ST*GCL*GTK ...
       + PZ*ST*GCL*GTNA + PY*PZ*GK*GCL);
    VB  = -(JNAK*FAR*ST*GTK + JNAK*FAR*ST*GTNA ...
       + PZ*JNAK*FAR*GCL - PY*PZ*GK*GCL*VY ...
       + PZ*ST*GCL*GTNA*VTX - PY*GK*ST*GTK*VY ...
       - PY*GK*ST*GTNA*VY + PZ*ST*GCL*GTK*VTY ...
       + PZ*ST*GCL*GTK*VZ + PZ*ST*GCL*GTNA*VZ)...
       /(PY*GK*ST*GTK + PY*GK*ST*GTNA + PZ*ST*GCL*GTK ...
       + PZ*ST*GCL*GTNA + PY*PZ*GK*GCL);
    VT  = VA - VB;
    JZ  = GCL * PZ * ( VA + VZ ) / FAR;
    JTX = GTNA * ST * ( VT - VTX ) /FAR;
    JTY = GTK * ST * ( VT - VTY ) / FAR;
    JY  = GK * PY * ( VB - VY ) / FAR;
% ---------------------------------------------------
% EQUATIONS OF THE SYSTEM
% ---------------------------------------------------

    F(1) = NU * JW;
    F(2) = LAMBDA_I * ( JNK1 - 3 * JNAK + J1 - J4 - JW * (XI * E_I) );
    F(3) = LAMBDA_I * ( JNK1 + 2 * JNAK - JY - JW * (YI * E_I) );
    F(4) = LAMBDA_I * ( 2 * JNK1 + J4 + JZ - JW * (ZI * E_I) );
    F(5) = LAMBDA_I * ( JB - 2 * J4 - JW * (BI * E_I) );
    F(6) = LAMBDA_L * ( JTX - QT * (XL * E_L));
    F(7) = LAMBDA_L * ( JTY - QT * (YL * E_L));
    F(8) = ((J_IPR + J_RYR)*(CE-CA)-J_SERCA+DELTA*(J_IN-J_PM))*(1/V)-NU*JW*CA/V; 
    F(9) = (IP_PLC-IP_DEG)*(1/V)- NU * JW*IP/V;
    F(10) = (H_INF-H)/TAU;
    F(11) = (G_INF-G)/TAU_RYR;
    F(12) = (DELTA*(J_IN-J_PM))*(1/V)-NU * JW*CT/V;
%     F(12) = SCALE_HY*LAMBDA_I * ( JB - J1 - JW * HI/SCALE_HY* E_I )/ V;
    
    F = F';
end