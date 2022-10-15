%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author: Ivan Braz Scarpa Pereira                                    %%%
%%% Responsible professor: Antonelli Monti, Dr.                         %%%
%%% Supervisors: Rodrigo Teixeira Pinto, Dr.                            %%%
%%%              Epameinondas Kontos, Dr.                               %%%
%%%              Asimenia Korompili, Msc.                               %%%
%%% This is an algorithm designed based on the code used at the Master  %%%
%%% Thesis "Optimization of Control Parameters of MMC converters in HVDC%%% 
%%% transmission applications" done in a collaboration work between RWTH%%%
%%% Aachen University and HVDC Control and Protection 1 from Siemens    %%%
%%% Energy, Erlangen, Germany                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
syms Theta_0 V_gd0 V_gq0 I_cd0g I_cd0 I_cq0 I_cq0g I_gd0 I_gq0 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 b_31 b_32 b_33
syms Ki_PLL Kp_PLL La Lt Ra Rt w0 VDCbase I_uq0 I_lq0 I_ud0 I_ld0 Req Leq V_cq0 V_cd0
syms Kp_Q Ki_Q Ki_DC Ki_P Kp_DC Kp_P I_g_sd0 I_g_sq0 V_gd0 V_gq0 Sbase VACMAG Ki_decd Ki_decq Kp_decd Kp_decq wLpu
syms Kid_CCSC Kiq_CCSC Kpd_CCSC Kpq_CCSC Vzd0 Vzq0 Vzd0_2c Vzq0_2c Leq Vtr1 Vtr2 Vsd0ref Vsq0ref Lg Rg VDCref

%% Reference transformation
% This section refers to the reference transformation modules as a function of the mismatch of the dq rotating
% frames between primary and secondary of the Y-Delta trafo
    VgToVcc    = [cos(Theta_0-pi()/6),sin(Theta_0-pi()/6),-V_gd0*sin(Theta_0-pi()/6)+V_gq0*cos(Theta_0-pi()/6);
        -sin(Theta_0-pi()/6),cos(Theta_0-pi()/6),-cos(Theta_0-pi()/6)*V_gd0-V_gq0*sin(Theta_0-pi()/6)];
    IccToIc2c = [cos(-3*Theta_0),sin(-3*Theta_0),3*sin(-3*Theta_0)*I_cd0g-3*cos(-3*Theta_0)*I_cq0g;
        -sin(-3*Theta_0),cos(-3*Theta_0),3*(cos(-3*Theta_0)*I_cd0g+sin(-3*Theta_0)*I_cq0g)];
    IsgToIsc   = [cos(Theta_0-pi()/6),sin(Theta_0-pi()/6),-I_gd0*sin(Theta_0-pi()/6)+I_gq0*cos(Theta_0-pi()/6);
        -sin(Theta_0-pi()/6),cos(Theta_0-pi()/6),-cos(Theta_0-pi()/6)*I_gd0-sin(Theta_0-pi()/6)*I_gq0];
    IscToIsg   = [cos(Theta_0-pi()/6),-sin(Theta_0-pi()/6),-I_cd0*sin(Theta_0-pi()/6)-I_cq0*cos(Theta_0-pi()/6);
        sin(Theta_0-pi()/6),cos(Theta_0-pi()/6),cos(Theta_0-pi()/6)*I_cd0-sin(Theta_0-pi()/6)*I_cq0];

    FeedbackVgdq = sqrt(3)* [Req,-wLpu,Leq,0,-Leq*I_cq0,1,0;
                    wLpu, Req,0,Leq,Leq*I_cd0,0,1];

    VcToVg= [cos(Theta_0-pi/6),-sin(Theta_0-pi/6),-sin(Theta_0-pi/6)*V_gd0-cos(Theta_0-pi/6)*V_gq0;
             sin(Theta_0-pi/6),cos(Theta_0-pi/6),cos(Theta_0-pi/6)*V_gd0-sin(Theta_0-pi/6)*V_gq0
    ];
    
    Vsd0refg = Vsd0ref;
    Vsq0refg = Vsq0ref;
D_grid = [-Rg+Lg*(Req+Rg)/(Leq+Lg),0,Lg/(Leq+Lg),0,(-1+sqrt(3)*Lg/(Leq+Lg)),0;
          0,-Rg+Lg*(Req+Rg)/(Leq+Lg),0,Lg/(Leq+Lg),0,(-1+sqrt(3)*Lg/(Leq+Lg));
];

VscPhaseToVscLine = [cos(Theta_0-pi/6),sin(Theta_0-pi/6),-sin(Theta_0-pi/6)*Vsd0refg+cos(Theta_0-pi/6)*Vsq0refg;
             -sin(Theta_0-pi/6),cos(Theta_0-pi/6),-cos(Theta_0-pi/6)*Vsd0refg-sin(Theta_0-pi/6)*Vsq0refg
    ];

VscToVsg = [cos(Theta_0-pi/6),sin(Theta_0-pi/6),-sin(Theta_0-pi/6)*Vsd0refg+cos(Theta_0-pi/6)*Vsq0refg;
             -sin(Theta_0-pi/6),cos(Theta_0-pi/6),-cos(Theta_0-pi/6)*Vsd0refg-sin(Theta_0-pi/6)*Vsq0refg
    ];

InterMatrixForIs_g = 1/sqrt(3)*[-Rg+Lg*(Req+Rg)/(Leq+Lg),0;
                    0,-Rg+Lg*(Req+Rg)/(Leq+Lg)
                    ]*IscToIsg;
InterMatrixForVs_g = sqrt(3)*[Lg/(Leq+Lg),0;0,Lg/(Leq+Lg)]*VscToVsg;

InterMatrixGeral_g = sqrt(3)*[InterMatrixForIs_g(1,1),InterMatrixForIs_g(1,2),InterMatrixForVs_g(1,1),InterMatrixForVs_g(1,2),D_grid(1,5),D_grid(1,6),InterMatrixForIs_g(1,3)+InterMatrixForVs_g(1,3);
                    InterMatrixForIs_g(2,1),InterMatrixForIs_g(2,2),InterMatrixForVs_g(2,1),InterMatrixForVs_g(2,2),D_grid(2,5),D_grid(2,6),InterMatrixForIs_g(2,3)+InterMatrixForVs_g(2,3);
];

InterMatrixGeral_c_semTheta = ([VgToVcc(1,1),VgToVcc(1,2);VgToVcc(2,1),VgToVcc(2,2)]*InterMatrixGeral_g);

InterMatrixGeral_c_semTheta(1,7) = InterMatrixGeral_c_semTheta(1,7)+ VgToVcc(1,3);
InterMatrixGeral_c_semTheta(2,7) = InterMatrixGeral_c_semTheta(2,7)+ VgToVcc(2,3);
InterMatrixGeral = sqrt(3)*InterMatrixGeral_c_semTheta;

%% Interconnection matrices from Component Connection Method

% For a better understading see Appendix C of the Master Thesis
% L11 interconnection matrix for Inverter
L11_1a =[InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0;...
                                              0,0,0,0,0,0,0,0,1,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,1,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                                              0,0,0,0,0,0,1,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0;
3*sin(-3*Theta_0)*I_cd0g-3*cos(-3*Theta_0)*I_cq0g,0,0,0,cos(-3*Theta_0),sin(-3*Theta_0),0,0,0,0,0,0,0,0;
3*(cos(-3*Theta_0)*I_cd0g+sin(-3*Theta_0)*I_cq0g),0,0,0,-sin(-3*Theta_0),cos(-3*Theta_0),0,0,0,0,0,0,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                                              1,0,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,1,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,1,0,0,0,0;
        ];
    %L12 interconnection matrix for Inverter
        L12_1a = [InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,1,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,1,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;

    ];
    %L12 interconnection matrix for Rectifier
L11_2a =[InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;...
                                              0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                                              0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
3*sin(-3*Theta_0)*I_cd0g-3*cos(-3*Theta_0)*I_cq0g,0,0,0,cos(-3*Theta_0),sin(-3*Theta_0),0,0,0,0,0,0,0,0,0;
3*(cos(-3*Theta_0)*I_cd0g+sin(-3*Theta_0)*I_cq0g),0,0,0,-sin(-3*Theta_0),cos(-3*Theta_0),0,0,0,0,0,0,0,0,0;
                                              0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
                                              1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
                                              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
         InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0,0;
         InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;
                                              0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
                                              0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        ];
    
    %L12 interconnection matrix for Rectifier
    L12_2a = [InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,1,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,1,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,0,0,0,0,0,0;
                  InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,1,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;
                          0,           0,0,0,0,0,0,0,0,0,0;

    ];

% Definition of P and Q as a function of the small signal parameters
% Their definition can be found at the master thesis "Optimization of MMC
% control parameters using MOOA"
P_theta = InterMatrixGeral(1,7)*I_cd0 + InterMatrixGeral(2,7)*I_cq0;
P_Isd   = InterMatrixGeral(1,1)*I_cd0 + InterMatrixGeral(2,1)*I_cq0 + V_gd0;
P_Isq   = InterMatrixGeral(1,2)*I_cd0 + InterMatrixGeral(2,2)*I_cq0 + V_gq0;
P_Vsd   = InterMatrixGeral(1,3)*I_cd0 + InterMatrixGeral(2,3)*I_cq0;
P_Vsq   = InterMatrixGeral(1,4)*I_cd0 + InterMatrixGeral(2,4)*I_cq0;
P_egd   = InterMatrixGeral(1,5)*I_cd0 + InterMatrixGeral(2,5)*I_cq0;
P_egq   = InterMatrixGeral(1,6)*I_cd0 + InterMatrixGeral(2,6)*I_cq0;

Q_theta = InterMatrixGeral(1,7)*I_cq0 + InterMatrixGeral(2,7)*I_cd0;
Q_Isd   = InterMatrixGeral(1,1)*I_cq0 + InterMatrixGeral(2,1)*I_cd0 + V_gq0;
Q_Isq   = InterMatrixGeral(1,2)*I_cq0 + InterMatrixGeral(2,2)*I_cd0 + V_gd0;
Q_Vsd   = InterMatrixGeral(1,3)*I_cq0 + InterMatrixGeral(2,3)*I_cd0;
Q_Vsq   = InterMatrixGeral(1,4)*I_cq0 + InterMatrixGeral(2,4)*I_cd0;
Q_egd   = InterMatrixGeral(1,5)*I_cq0 + InterMatrixGeral(2,5)*I_cd0;
Q_egq   = InterMatrixGeral(1,6)*I_cq0 + InterMatrixGeral(2,6)*I_cd0;

%L21 interconnection matrix for Inverter
L_21 = [0,0,1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0;
InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0;
InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0
        VscPhaseToVscLine(1,3),0,0,0,0,0,0,0,0,0,0,0,VscPhaseToVscLine(1,1),VscPhaseToVscLine(1,2);
        -VscPhaseToVscLine(2,3),0,0,0,0,0,0,0,0,0,0,0,-VscPhaseToVscLine(2,1),-VscPhaseToVscLine(2,2);
        P_theta,0,P_Isd,P_Isq,0,0,0,0,P_Vsd,P_Vsq,0,0,0,0;
        Q_theta,0,Q_Isd,Q_Isq,0,0,0,0,Q_Vsd,Q_Vsq,0,0,0,0;
        0,0,0,0,0,0,1,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0,0,0,0;
];

%L21 interconnection matrix for Rectifier
L_21b = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;
InterMatrixGeral(1,7),0,InterMatrixGeral(1,1),InterMatrixGeral(1,2),0,0,0,0,InterMatrixGeral(1,3),InterMatrixGeral(1,4),0,0,0,0,0;
InterMatrixGeral(2,7),0,InterMatrixGeral(2,1),InterMatrixGeral(2,2),0,0,0,0,InterMatrixGeral(2,3),InterMatrixGeral(2,4),0,0,0,0,0;
        VscPhaseToVscLine(1,3),0,0,0,0,0,0,0,0,0,0,0,VscPhaseToVscLine(1,1),VscPhaseToVscLine(1,2),0;
        -VscPhaseToVscLine(2,3),0,0,0,0,0,0,0,0,0,0,0,-VscPhaseToVscLine(2,1),-VscPhaseToVscLine(2,2),0;
];

%L22 interconnection matrix for Inverter
L_22 = [0,0,0,0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,0,0,0,1;
 InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,1,0,0,0,0,0;
 InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,0,1,0,0;
        P_egd,P_egq,0,0,0,0,0,0,0,0,0;
        Q_egd,Q_egq,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,0;
];

%L22 interconnection matrix for Rectifier
L_22DC = [0,0,0,0,0,0,0,0,0,1,0;
        0,0,0,0,0,0,0,0,0,0,1;
 InterMatrixGeral(1,5),InterMatrixGeral(1,6),0,0,0,1,0,0,0,0,0;
 InterMatrixGeral(2,5),InterMatrixGeral(2,6),0,0,0,0,1,0,0,0,0;
        0,0,0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,0,1,0,0;
];

%% Modelling of the Control and Converter modules
%Calculation of the DC voltage (DC grid)

A_Vdc = [0];
B_Vdc  = 1/VDCref*[1, I_gd0, I_gq0, V_gd0, V_gq0];
C_Vdc  = [1];
D_Vdc = [0,0,0,0,0];

% PLL
A_PLL = [0,0;Ki_PLL,0];
B_PLL = [-1;-Kp_PLL];
C_PLL = [0,1;Ki_PLL,0];
D_PLL = [0;-Kp_PLL];
    
%Converter Model

A_conv = [-Req/Leq,w0,0,0;
          -w0,-Req/Leq,0,0;
          0,0,-Ra/La,w0;
          0,0,-w0,-Ra/La
          ];
B_conv = [-1/Leq,0,1/(Leq),0,0,0,I_cq0;
          0,-1/Leq,0,1/(Leq),0,0,-I_cd0;
          0,0,0,0,-1/La,0,I_cq0g;
          0,0,0,0,0,-1/La,-I_cd0g
];

C_conv = [1,0,0,0;
          0,1,0,0;
          0,0,1,0;
          0,0,0,1
];
       
D_conv = [0,0,0,0,0,0,0;
          0,0,0,0,0,0,0;
          0,0,0,0,0,0,0;
          0,0,0,0,0,0,0;
];
    
    
%Real Power loop

A_P = [0];
B_P = [1, -I_gd0, -I_gq0, -V_gd0, -V_gq0];
C_P = 1/VACMAG*[Ki_P];
D_P = 1/VACMAG*[Kp_P, -Kp_P*I_gd0,-Kp_P*I_gq0, -Kp_P*V_gd0, -Kp_P*V_gq0];

%DC Voltage Loop

A_DC = [0];
B_DC = [1, -1];
C_DC = [Ki_DC];
D_DC = [Kp_DC, -Kp_DC];

%Reactive Power Loop

A_Q = [0];
B_Q = [1, -I_gq0, -I_gd0, -V_gq0, -V_gd0];
C_Q = [Ki_Q/VACMAG];
D_Q = [Kp_Q/(VACMAG), -Kp_Q*I_gq0/VACMAG,-Kp_Q*I_gd0/VACMAG,-Kp_Q*V_gq0/VACMAG, -Kp_Q*V_gd0/VACMAG];

%Active and Reactive Together

A_PQ = [0,0;0,0];
B_PQ = [1,0, -I_gd0, -I_gq0, -V_gd0, -V_gq0;
        0,1, -I_gq0, -I_gd0, -V_gq0, -V_gd0
];
C_PQ = 1/VACMAG*[Ki_P,0;0,Ki_Q];
D_PQ = 1/VACMAG*[Kp_P,0, -Kp_P*I_gd0,-Kp_P*I_gq0, -Kp_P*V_gd0, -Kp_P*V_gq0;
                    0,Kp_Q, -Kp_Q*I_gq0,-Kp_Q*I_gd0,-Kp_Q*V_gq0, -Kp_Q*V_gd0
];

% Phase Current Controller

A_dec = [-Ki_decd/Kp_decd,0;0,-Ki_decq/Kp_decq];
B_dec = [1/(Kp_decd),0,-1/Kp_decd,0,0,0,0,wLpu/Kp_decd,Leq*I_cq0/Kp_decd;0,1/(Kp_decq),0,-1/(Kp_decq),0,0,-wLpu/Kp_decq,0,-Leq*I_cd0/Kp_decq];
C_dec = [-Ki_decd,0;0,-Ki_decq];
D_dec = [1,0,0,0,-Kp_decd,0,Kp_decd,wLpu,Leq*I_cq0;0,1,0,0,0,-Kp_decq,-wLpu,Kp_decq,-Leq*I_cd0];

a= cos(3*Theta_0);
b= sin(3*Theta_0);
c= -3*sin(3*Theta_0)*Vzd0_2c+3*cos(3*Theta_0)*Vzq0_2c;
d= -sin(3*Theta_0);
e= cos(3*Theta_0);
f= -3*cos(3*Theta_0)*Vzd0_2c -3*sin(3*Theta_0)*Vzq0_2c;

%CCSC
A_CCSC = [0,0;0,0];
B_CCSC = [-1,0,0,0,0,0;0,-1,0,0,0,0];
C_CCSC = [
          -a*Kid_CCSC,-b*Kiq_CCSC;
          -d*Kid_CCSC,-e*Kiq_CCSC;
          0,0;
          0,0
           ];
       
 
D_CCSC = [a*Kpd_CCSC-2*w0*La*b,2*w0*La*a+b*Kpq_CCSC,2*La*I_cq0g*a-2*La*I_cd0g*b,c,0,0;
          d*Kpd_CCSC-2*w0*La*e,2*d*w0*La+Kpq_CCSC*e,2*d*La*I_cq0g-2*La*I_cd0g*e,f,0,0;
          0                   ,0                   ,0                          ,0,1,0;
          0                   ,0                   ,0                          ,0,0,1
];

%% Overall state space matrices for Rectifier and Inverter
    % A_star definition for Inverter
    A_star= sym(zeros(size(A_PLL,1)+size(A_CCSC,1)+size(A_dec,1)+size(A_Q,1)+size(A_P,1)+size(A_conv,1),size(A_PLL,2)+size(A_CCSC,2)+size(A_dec,2)+size(A_Q,2)+size(A_P,2)+size(A_conv,2)));
    A_star(1:2,1:2)=A_PLL;
    A_star(3:6,3:6)=A_conv;
    A_star(7,7)= A_P;
    A_star(8,8)= A_Q;
    A_star(9:10,9:10)=A_dec;
    A_star(11:12,11:12)=A_CCSC;

    %B_star definition for Inverter
    B_star= sym(zeros(size(B_PLL,1)+size(B_CCSC,1)+size(B_dec,1)+size(B_Q,1)+size(B_P,1)+size(B_conv,1),size(B_PLL,2)+size(B_CCSC,2)+size(B_dec,2)+size(B_Q,2)+size(B_P,2)+size(B_conv,2)));
    B_star(1:2,1)=B_PLL;
    B_star(3:6,2:8)=B_conv;
    B_star(7,9:13)= B_P;
    B_star(8,14:18)= B_Q;
    B_star(9:10,19:27)=B_dec;
    B_star(11:12,28:33)=B_CCSC;
    % 
    %C_star definition for Inverter
    C_star= sym(zeros(size(C_PLL,1)+size(C_CCSC,1)+size(C_dec,1)+size(C_Q,1)+size(C_P,1)+size(C_conv,1),size(C_PLL,2)+size(C_CCSC,2)+size(C_dec,2)+size(C_Q,2)+size(C_P,2)+size(C_conv,2)));
    C_star(1:2,1:2)=C_PLL;
    C_star(3:6,3:6)=C_conv;
    C_star(7,7)= C_P;
    C_star(8,8)= C_Q;
    C_star(9:10,9:10)=C_dec;
    C_star(11:14,11:12)=C_CCSC;
    % 
    %D_star definition for Inverter
    D_star= sym(zeros(size(D_PLL,1)+size(D_CCSC,1)+size(D_dec,1)+size(D_Q,1)+size(D_P,1)+size(D_conv,1),size(D_PLL,2)+size(D_CCSC,2)+size(D_dec,2)+size(D_Q,2)+size(D_P,2)+size(D_conv,2)));
    D_star(1:2,1)=D_PLL;
    D_star(3:6,2:8)=D_conv;
    D_star(7,9:13)= D_P;
    D_star(8,14:18)= D_Q;
    D_star(9:10,19:27)=D_dec;
    D_star(11:14,28:33)=D_CCSC;

    % A_star definition for Rectifier
    A_star_V= sym(zeros(size(A_PLL,1)+size(A_CCSC,1)+size(A_dec,1)+size(A_Q,1)+size(A_DC,1)+size(A_conv,1),size(A_PLL,2)+size(A_CCSC,2)+size(A_dec,2)+size(A_Q,2)+size(A_DC,2)+size(A_conv,2)));
    A_star_V(1:2,1:2)=A_PLL;
    A_star_V(3:6,3:6)=A_conv;
    A_star_V(7,7)= A_DC;
    A_star_V(8,8)= A_Q;
    A_star_V(9:10,9:10)=A_dec;
    A_star_V(11:12,11:12)=A_CCSC;
    A_star_V(13,13)=A_Vdc;

    %B_star definition for Rectifier
    B_star_V= sym(zeros(size(B_PLL,1)+size(B_CCSC,1)+size(B_dec,1)+size(B_Q,1)+size(B_DC,1)+size(B_conv,1),size(B_PLL,2)+size(B_CCSC,2)+size(B_dec,2)+size(B_Q,2)+size(B_DC,2)+size(B_conv,2)));
    B_star_V(1:2,1)=B_PLL;
    B_star_V(3:6,2:8)=B_conv;
    B_star_V(7,9:10)= B_DC;
    B_star_V(8,11:15)= B_Q;
    B_star_V(9:10,16:24)=B_dec;
    B_star_V(11:12,25:30)=B_CCSC;
    B_star_V(13,31:35)=B_Vdc;
    % 
    %C_star definition for Rectifier
    C_star_V= sym(zeros(size(C_PLL,1)+size(C_CCSC,1)+size(C_dec,1)+size(C_Q,1)+size(C_DC,1)+size(C_conv,1),size(C_PLL,2)+size(C_CCSC,2)+size(C_dec,2)+size(C_Q,2)+size(C_DC,2)+size(C_conv,2)));
    C_star_V(1:2,1:2)=C_PLL;
    C_star_V(3:6,3:6)=C_conv;
    C_star_V(7,7)= C_DC;
    C_star_V(8,8)= C_Q;
    C_star_V(9:10,9:10)=C_dec;
    C_star_V(11:14,11:12)=C_CCSC;
    C_star_V(15,13)=C_Vdc;
    % 
    %D_star definition for Rectifier
    D_star_V= sym(zeros(size(D_PLL,1)+size(D_CCSC,1)+size(D_dec,1)+size(D_Q,1)+size(D_DC,1)+size(D_conv,1),size(D_PLL,2)+size(D_CCSC,2)+size(D_dec,2)+size(D_Q,2)+size(D_DC,2)+size(D_conv,2)));
    D_star_V(1:2,1)=D_PLL;
    D_star_V(3:6,2:8)=D_conv;
    D_star_V(7,9:10)= D_DC;
    D_star_V(8,11:15)= D_Q;
    D_star_V(9:10,16:24)=D_dec;
    D_star_V(11:14,25:30)=D_CCSC;
    D_star_V(15,31:35)=D_Vdc;

    %State Space overall definition

    %P control
    Atotal_Pa = sym(A_star + B_star*L11_1a*inv(eye(14)-D_star*L11_1a)*C_star);
    Btotal_Pa = sym(B_star*L11_1a*inv(eye(14)-D_star*L11_1a)*D_star*L12_1a+B_star*L12_1a);
    Ctotal_Pa = sym(L_21*inv(eye(14)-D_star*L11_1a)*C_star);
    Dtotal_Pa = sym(L_21*inv(eye(14)-D_star*L11_1a)*D_star*L12_1a+L_22);


    %VDC control
    Atotal_DCa = sym(A_star_V + B_star_V*L11_2a*inv(eye(15)-D_star_V*L11_2a)*C_star_V);
    Btotal_DCa = sym(B_star_V*L11_2a*inv(eye(15)-D_star_V*L11_2a)*D_star_V*L12_2a+B_star_V*L12_2a);
    Ctotal_DCa = sym(L_21b*inv(eye(15)-D_star_V*L11_2a)*C_star_V);
    Dtotal_DCa = sym(L_21b*inv(eye(15)-D_star_V*L11_2a)*D_star_V*L12_2a+L_22DC);