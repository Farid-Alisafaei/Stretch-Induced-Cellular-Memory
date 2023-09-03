%
% Matrix Softening Controls Stretch-Induced Cellular Memory and Fibroblast Activation
%
% For any questions regarding this file, please contact Farid Alisafaei (farid.alisafaei@njit.edu)
%
%--------------------------------------------------------------------------
%
clear
clc
%
% time (hr)
tspan = [0 24] ;
%
% initial guess for epsilonM (matrix strain)
% matrix strain is denoted as ε_m in the manuscript
y0 = 0 ;
%
% stiffness of ECM in compression (kPa)
% denoted as ε_c in the manuscript
% refer to "SI section 1.2"
E_MT = 9 ;
%
% cell initial contractility (kPa)
% denoted as σ_c^0 in the manuscript
% see table S1
rho_0 = 6.2  ;
%
% cell chemo-mechanical feedback parameter 1/(kPa)
% denoted as α in the manuscript
% see table S1 
alpha = 0.45 ;
%
% cell chemical stiffness parameter 1/(kPa)
% denoted as β in the manuscript
% see table S1 
beta  = 0.5 ;
%
% viscosity of the matrix in tension (kPa/hr) 
% denoted as η in the manuscript
% refer to "SI section 1.1"
yita  = 170 ;
%
% stablity criterion
% described in "SI section 1.5"
if ((beta-alpha) < 0)
    disp('ERROR')
    return
elseif ((beta-alpha) == 0)
    disp('ERROR')
    return
elseif (E_MT*alpha-1 < 0)
    disp('ERROR')
    return
elseif (E_MT*alpha-1 == 0)
    disp('ERROR')
    return
else
    disp('stablity criterion satisfied')
end
disp('----------------------------')
%
% defining the constitutive equation for the matrix element in tension
% refer to "SI section 1.1" and equation S2
epsilon1        = 0.200                            ; % denoted as ε_1 in equation S2
epsilon2        = 0.670                            ; % denoted as ε_2 in equation S2
transitionPoint = epsilon2 - epsilon1              ; % denoted as ε_t in equation S2
yita1           = pi/epsilon2                      ; % denoted as η_1 in equation S2
ksi1            = yita1*cot(yita1*epsilon1)        ; % denoted as ξ_1 in equation S2 
epsilon4        = 1.10                             ; % denoted as ε_4 in equation S2
yita2           = pi / epsilon4                    ; % denoted as ξ_2 in equation S2 
ksi2            = yita2*cot(yita2*transitionPoint) ; % denoted as ξ_2 in equation S2
YY              = []                               ;
%
% maximum tissue strain
epsilon_tot_max = 0.35 ;
%
% tissue strain
% denoted as ε in the manuscript
Increament  = 0.01                           ;
epsilon_tot = [0:Increament:epsilon_tot_max] ;
%
% length of the tissue (mm)
L = 13.3 ;
%
% width of the tissue (mm)
A = 1.0 ;
%
% thickness of the tissue (mm)
B = 0.1 ;
%
% cross section area of the tissue (mm^2)
Area = A * B ;
%
%--------------------------------------------------------------------------
%
for i = 1:length(epsilon_tot)
    %
    % solving for the strain of the matrix element in tension (denoted as ε_m in the manuscript)
    [ t , epsilon ] = ode45(@(t,epsilon) odefun( t , epsilon , E_MT , alpha , beta , epsilon_tot(i) , rho_0 , yita, yita1 , yita2 , ksi1 , ksi2 , epsilon2 , transitionPoint ), tspan , y0 );
    %
    Index(i) = length(epsilon) ;
    %
    % time (hr)
    Time( 1:Index(i) , i ) = t ;
    %
    % strain of the matrix element in tension 
    % denoted as ε_m in the manuscript
    % the tensor rows display ε_m at different tissue strains e
    % the tensor columns display ε_m at different time points t
    epsilon_M( 1:Index(i) , i ) = epsilon ;
    %
    % cell strain 
    % denoted as ε_c in the manuscript
    % see equation S10
    % the tensor rows display ε_c at different tissue strains ε
    % the tensor columns display ε_c at different time points t
    epsilon_C( 1:Index(i) , i ) = epsilon_tot(i) - epsilon_M( 1:Index(i) , i ) ;
    %
    % cell contractility (kPa)
    % denoted as σ_c in the manuscript
    % see equation S8
    % the tensor rows display σ_c at different tissue strains ε
    % the tensor columns display σ_c at different time points t
    rho( 1:Index(i) , i )       = ( (E_MT*alpha- 1)/(beta-alpha) ) .* epsilon_C( 1:Index(i) , i )   +   beta*rho_0/(beta-alpha) ;
    %
    % stress (kPa)
    % denoted as σ in the manuscript
    % see equation S7
    % the tensor rows display σ at different tissue strains ε
    % the tensor columns display σ at different time points t
    sigma( 1:Index(i) , i )     = rho( 1:Index(i) , i ) + E_MT * epsilon_C( 1:Index(i) , i ) ;
    %
    % non-viscous stress of the matrix element in tension (kPa)
    % denoted as σ_s in the manuscript
    % see equation S2
    % the tensor rows display σ_c at different tissue strains ε
    % the tensor columns display σ_c at different time points t
    sigma_S( 1: Index(i) , i)   = ((40*exp(ksi1/yita1*atan(yita1/ksi1))*sqrt(1+(ksi1/yita1)^2)*exp(-ksi1*(epsilon2-epsilon_M( 1:Index(i) , i ))).*sin(yita1*(epsilon2-epsilon_M( 1:Index(i) , i )))).*(epsilon_M( 1:Index(i) , i )<=transitionPoint)+(40*exp(ksi2/yita2*atan(yita2/ksi2))*sqrt(1+(ksi2/yita2)^2)*exp(-ksi2*(epsilon_M( 1:Index(i) , i ))).*sin(yita2*(epsilon_M( 1:Index(i) , i )))).*(epsilon_M( 1:Index(i) , i )>transitionPoint));
    %
    % viscous stress of the matrix element in tension (kPa)
    % denoted as σ_η in the manuscript
    % see equation S3
    % the tensor rows display σ_η at different tissue strains ε
    % the tensor columns display σ_η at different time points t
    sigma_yita( 1:Index(i) , i ) = sigma( 1:Index(i) , i ) - sigma_S( 1:Index(i) , i ) ;
    %
    % tissue strain is 0%
    if ( epsilon_tot(i) == 0 )
        Index_00 = i ;
    end
    %
    % tissue strain is 5%
    if ( epsilon_tot(i) == 0.05 )
        Index_05 = i ;
    end
    %
    % tissue strain is 15%
    if ( epsilon_tot(i) == 0.15 )
        Index_15 = i ;
    end
    %
    % tissue strain is 30%
    if ( epsilon_tot(i) == 0.30 )
        Index_30 = i ;
    end
    % 
end
%
%--------------------------------------------------------------------------
%
% constitutive model for the matrix in tension
%
% ε_m
epsilon_MM = linspace(0,epsilon4,100) ;
%
% σ_s as a function of ε_m
% see equation S2
sigma_SS = ((40*exp(ksi1/yita1*atan(yita1/ksi1))*sqrt(1+(ksi1/yita1)^2)*exp(-ksi1*(epsilon2-epsilon_MM))         .*sin(yita1*(epsilon2-epsilon_MM)))         .*(epsilon_MM<=transitionPoint)        +(40*exp(ksi2/yita2*atan(yita2/ksi2))*sqrt(1+(ksi2/yita2)^2)*exp(-ksi2*(epsilon_MM))            .*sin(yita2*(epsilon_MM)))            .*(epsilon_MM>transitionPoint));
%
%--------------------------------------------------------------------------
%
% tissue stress σ
%
figure
plot( Time( 1:Index(Index_00) , Index_00 ) , sigma( 1:Index(Index_00) , Index_00 ) , 'Color' , [0      0      1     ] )
hold on
plot( Time( 1:Index(Index_15) , Index_15 ) , sigma( 1:Index(Index_15) , Index_15 ) , 'Color' , [1      0      0     ] )
hold on
plot( Time( 1:Index(Index_30) , Index_30 ) , sigma( 1:Index(Index_30) , Index_30 ) , 'Color' , [0.2549 0.5804 0.1294] )
xlim([0, 24])
xlabel('Time (h)')
ylabel('Tissu stress, σ (kPa)')
%
for i = 1:length(epsilon_tot)
    X_epsilon_tot(i) = epsilon_tot(i)        ;
    Y_sigma(i)       = sigma( Index(i) , i ) ;
end
figure
plot( X_epsilon_tot , Y_sigma , 'k' )
xlim([0, 0.3])
xlabel('Tissue strain, ε')
ylabel('Tissue stress, σ (kPa)')
%
%--------------------------------------------------------------------------
%
% cell stress σ_c
%
figure
plot( Time( 1:Index(Index_00) , Index_00 ) , rho( 1:Index(Index_00) , Index_00 ) , 'Color' , [0      0      1     ] )
hold on
plot( Time( 1:Index(Index_15) , Index_15 ) , rho( 1:Index(Index_15) , Index_15 ) , 'Color' , [1      0      0     ] )
hold on
plot( Time( 1:Index(Index_30) , Index_30 ) , rho( 1:Index(Index_30) , Index_30 ) , 'Color' , [0.2549 0.5804 0.1294] )
xlim([0, 24])
xlabel('Time (h)')
ylabel('Cell stress, σ_c (kPa)')
%
for i = 1:length(epsilon_tot)
    Y_rho(i)         = rho( Index(i) , i ) ;
end
figure
plot( X_epsilon_tot , Y_rho , 'k' )
xlim([0, 0.3])
xlabel('Tissue strain, ε')
ylabel('Cell stress, σ_c (kPa)')
%
%--------------------------------------------------------------------------
%
% cell strain ε_c
%
figure
plot( Time( 1:Index(Index_00) , Index_00 ) , epsilon_C( 1:Index(Index_00) , Index_00 ) , 'Color' , [0      0      1     ] )
hold on
plot( Time( 1:Index(Index_15) , Index_15 ) , epsilon_C( 1:Index(Index_15) , Index_15 ) , 'Color' , [1      0      0     ] )
hold on
plot( Time( 1:Index(Index_30) , Index_30 ) , epsilon_C( 1:Index(Index_30) , Index_30 ) , 'Color' , [0.2549 0.5804 0.1294] )
xlim([0, 24])
xlabel('Time (h)')
ylabel('cell strain, ε_c')
%
for i = 1:length(epsilon_tot)
    Y_epsilon_C(i) = epsilon_C( Index(i) , i ) ;
end
figure
plot( X_epsilon_tot ,Y_epsilon_C , 'k' )
xlim([0, 0.3])
xlabel('Tissue strain, ε')
ylabel('Cell strain, ε_c')
%
%--------------------------------------------------------------------------
%
% matrix strain ε_m
%
figure
plot( Time( 1:Index(Index_00) , Index_00 ) , epsilon_M( 1:Index(Index_00) , Index_00 ) , 'Color' , [0      0      1     ] )
hold on
plot( Time( 1:Index(Index_15) , Index_15 ) , epsilon_M( 1:Index(Index_15) , Index_15 ) , 'Color' , [1      0      0     ] )
hold on
plot( Time( 1:Index(Index_30) , Index_30 ) , epsilon_M( 1:Index(Index_30) , Index_30 ) , 'Color' , [0.2549 0.5804 0.1294] )
xlim([0, 24])
xlabel('Time (h)')
ylabel('Matrix strain, ε_m')
%
for i = 1:length(epsilon_tot)
    Y_epsilon_M(i) = epsilon_M( Index(i) , i ) ;
end
figure
plot( X_epsilon_tot , Y_epsilon_M , 'k' )
xlim([0, 0.3])
xlabel('Tissue strain, ε')
ylabel('Matrix strain, ε_m')
%
%--------------------------------------------------------------------------
%
% local matrix stiffness after 24h tissue stretching is 
% determined by calculating the slope of the stress-strain curve for the matrix element in tension
% see Figure 4B (right panel)
% denoted as E_tangent in the manuscript
m_end_00 = ( sigma_S( Index(Index_00+1) , Index_00+1 ) - sigma_S( Index(Index_00) , Index_00 ) ) / ( epsilon_M( Index(Index_00+1) , Index_00+1 ) - epsilon_M( Index(Index_00) , Index_00 ) ) ; 
m_end_15 = ( sigma_S( Index(Index_15+1) , Index_15+1 ) - sigma_S( Index(Index_15) , Index_15 ) ) / ( epsilon_M( Index(Index_15+1) , Index_15+1 ) - epsilon_M( Index(Index_15) , Index_15 ) ) ; 
m_end_30 = ( sigma_S( Index(Index_30+1) , Index_30+1 ) - sigma_S( Index(Index_30) , Index_30 ) ) / ( epsilon_M( Index(Index_30+1) , Index_30+1 ) - epsilon_M( Index(Index_30) , Index_30 ) ) ; 
%
disp(['Slope_24h_00% = ',num2str(m_end_00),' kPa'])
disp(['Slope_24h_15% = ',num2str(m_end_15),' kPa'])
disp(['Slope_24h_30% = ',num2str(m_end_30),' kPa'])
%
figure
plot( epsilon_M( Index(Index_00) , Index_00 ) , sigma_S( Index(Index_00) , Index_00 ) , 'o' , 'Color' , [0      0      1     ]  ) 
hold on
plot( epsilon_M( Index(Index_15) , Index_15 ) , sigma_S( Index(Index_15) , Index_15 ) , 'o' , 'Color' , [1      0      0     ]  ) 
hold on
plot( epsilon_M( Index(Index_30) , Index_30 ) , sigma_S( Index(Index_30) , Index_30 ) , 'o' , 'Color' , [0.2549 0.5804 0.1294]  ) 
hold on
plot( epsilon_MM , sigma_SS , 'k:' )
xlabel('Matrix strain, ε_m')
ylabel('Matrix stress (kPa)')
%
disp('----------------------------')
%
% 1h
%
for i = 1:Index(Index_00)
    Index_Time_00 = i ;
    if Time( i , Index_00 ) > 1
        break
    end
end
%
for i = 1:Index(Index_15)
    Index_Time_15 = i ;
    if Time( i , Index_15 ) > 1
        break
    end
end
%
for i = 1:Index(Index_30)
    Index_Time_30 = i ;
    if Time( i , Index_30 ) > 1
        break
    end
end
%
% local matrix stiffness after 1h stretching the tissue
% determined by calculating the slope of the stress-strain curve for the matrix element in tension
% see Figure 4B (left panel)
% denoted as E_tangent in the manuscript
m_beginning_00 = ( sigma_S( Index_Time_00+1 , Index_00+1 ) - sigma_S( Index_Time_00 , Index_00 ) ) / ( epsilon_M( Index_Time_00+1 , Index_00+1 ) - epsilon_M( Index_Time_00 , Index_00 ) ) ; 
m_beginning_15 = ( sigma_S( Index_Time_15+1 , Index_15+1 ) - sigma_S( Index_Time_15 , Index_15 ) ) / ( epsilon_M( Index_Time_15+1 , Index_15+1 ) - epsilon_M( Index_Time_15 , Index_15 ) ) ; 
m_beginning_30 = ( sigma_S( Index_Time_30+1 , Index_30+1 ) - sigma_S( Index_Time_30 , Index_30 ) ) / ( epsilon_M( Index_Time_30+1 , Index_30+1 ) - epsilon_M( Index_Time_30 , Index_30 ) ) ; 
%
disp(['Slope_1h_00% = ',num2str(m_beginning_00),' kPa'])
disp(['Slope_1h_15% = ',num2str(m_beginning_15),' kPa'])
disp(['Slope_1h_30% = ',num2str(m_beginning_30),' kPa'])
%
disp('----------------------------')
%
% tissues pre-strained for 1h
% see Figure 4B (left panel)
figure
plot( epsilon_M( Index_Time_00 , Index_00 ) , sigma_S( Index_Time_00        , Index_00 ) , 'o' , 'Color' , [0      0      1     ]  ) 
hold on
plot( epsilon_M( Index_Time_15 , Index_15 ) , sigma_S( Index_Time_15        , Index_15 ) , 'o' , 'Color' , [1      0      0     ]   ) 
hold on
plot( epsilon_M( Index_Time_30 , Index_30 ) , sigma_S( Index_Time_30        , Index_30 ) , 'o' , 'Color' , [0.2549 0.5804 0.1294]   ) 
hold on
plot( epsilon_MM , sigma_SS , 'k:' )
xlabel('Matrix strain, ε_m')
ylabel('Matrix stress (kPa)')
%
%--------------------------------------------------------------------------
%
% 1 h
Index_Time_1h = zeros(size(epsilon_tot)) ;
for i = 1:length(epsilon_tot)
    for j = 1:Index(Index_00)-1
        if ( (Time(j+1,i)>1)&&(Time(j,i)<1) )
            Index_Time_1h(i) = j+1 ;
        end
    end
end
%
% 1h
for i = 1:length(epsilon_tot)
    X_epsilon_tot(i)    = epsilon_tot(i)                    ;
    Y_sigma_1h(i)       = sigma(     Index_Time_1h(1) , i ) ;
    Y_rho_1h(i)         = rho(       Index_Time_1h(1) , i ) ;
    Y_epsilon_C_1h(i)   = epsilon_C( Index_Time_1h(1) , i ) ;
    Y_epsilon_M_1h(i)   = epsilon_M( Index_Time_1h(1) , i ) ;
end
%
% tissue stress (σ) vs tissue strain (ε) after 1h stretching
% see Figure 4C
figure
plot( X_epsilon_tot , Y_sigma_1h , 'k' )
xlabel('Tissue strain, ε')
ylabel('Tissu stress, σ (kPa)')
%
% cell stress (σ_c) vs tissue strain (ε) after 1h stretching
% see Figure 4C
figure
plot( X_epsilon_tot , Y_rho_1h , 'k' )
xlabel('Tissue strain, ε')
ylabel('Cell stress, σ_c (kPa)')
%
% cell strain (ε_c) vs tissue strain (ε) after 1h stretching
% see Figure 4C
figure
plot( X_epsilon_tot , Y_epsilon_C_1h , 'k' )
xlabel('Tissue strain, ε')
ylabel('Cell strain, ε_c')
%
% matrix strain (ε_m) vs tissue strain (ε) after 1h stretching
% see Figure 4C
figure
plot( X_epsilon_tot , Y_epsilon_M_1h , 'k' )
xlabel('Tissue strain, ε')
ylabel('Matrix strain, ε_m')
%
%--------------------------------------------------------------------------
%
% no matrix damage in the tissues stretched for 1h
% undamaged matrix stiffness (kPa)
% converting stiffness from kPa to N/m for the tissues stretched for 1h
E_M_Undamaged = m_beginning_15 ;
%
disp(['Elastic Modulus_1h_00% = ',num2str(m_beginning_00),' kPa  = ',num2str(m_beginning_00*Area/L),' N/m'])
disp(['Elastic Modulus_1h_15% = ',num2str(m_beginning_15),' kPa  = ',num2str(m_beginning_15*Area/L),' N/m'])
disp(['Elastic Modulus_1h_30% = ',num2str(m_beginning_30),' kPa  = ',num2str(m_beginning_30*Area/L),' N/m'])
%
% determine matrix stiffness for the tissues stretched for 24h 
% damage parameters
epsilon_M_1 = epsilon_M( Index(Index_00) , Index_00 ) ;  
epsilon_M_2 = epsilon4                                ;
D_1         = 1 - ( m_end_00 / E_M_Undamaged )        ;
D_2         = 1                                       ;
%
% damaged matrix stiffness for 24h loading with 0% strain (kPa)
D_00 = DamageFunction( epsilon_M_1 , epsilon_M_2 , D_1 , D_2 , epsilon_M( Index(Index_00) , Index_00 ) ) ;
E_M_00 = ( 1 - D_00 ) * E_M_Undamaged ; 
disp(['Elastic Modulus_24h_00% = ',num2str(E_M_00),' kPa  = ',num2str(E_M_00*Area/L),' N/m'])
%
% damaged matrix stiffness for 24h loading with 5% strain (kPa)
D_15 = DamageFunction( epsilon_M_1 , epsilon_M_2 , D_1 , D_2 , epsilon_M( Index(Index_15) , Index_15 ) ) ;
E_M_15 = ( 1 - D_15 ) * E_M_Undamaged ;
disp(['Elastic Modulus_24h_15% = ',num2str(E_M_15),' kPa  = ',num2str(E_M_15*Area/L),' N/m'])
%
% damaged matrix stiffness for 24h loading with 30% strain (kPa)
D_30 = DamageFunction( epsilon_M_1 , epsilon_M_2 , D_1 , D_2 , epsilon_M(   Index(Index_30) , Index_30 ) ) ;
E_M_30 = ( 1 - D_30 ) * E_M_Undamaged ;
disp(['Elastic Modulus_24h_30% = ',num2str(E_M_30),' kPa  = ',num2str(E_M_30*Area/L),' N/m'])
disp('----------------------------')
%
% matrix stiffness
% see Figure 2F
figure
Y_Stiffness = [E_M_00 E_M_15 E_M_30; m_beginning_00 m_beginning_15 m_beginning_30] * Area/L ;
bar(Y_Stiffness)
title('Matrix stiffness (N/m)')
%
%--------------------------------------------------------------------------
%
% baseline contractility (kPa)
% the minimum cell stress is obtained when the tissue is unconstrained and free to contract without any resistance (σ=0)
% refer to "SI section 1.3"
rho_min = E_MT * beta * rho_0 / ( E_MT * beta - 1 ) ; 
%
% cell contractiltiy at the end of 24h pre-straining (kPa) 
rho_prim_00_24h = mean( rho( Index(Index_00) , Index_00 ) ) ;  
rho_prim_15_24h = mean( rho( Index(Index_15) , Index_15 ) ) ;   
rho_prim_30_24h = mean( rho( Index(Index_30) , Index_30 ) ) ; 
%
% cell contractility at the end of 1h pre-straining (kPa)
rho_prim_00_1h  = mean( rho( Index_Time_00   , Index_00 ) ) ; 
rho_prim_15_1h  = mean( rho( Index_Time_15   , Index_15 ) ) ;
rho_prim_30_1h  = mean( rho( Index_Time_30   , Index_30 ) ) ;
%
% initial guess
y0 = [0 0 0 0] ;
%
% cell contractility after priming (Pa): stress that cells experience after removing the extrenal strain (Pa) 
rho_AfterPriming = rho_min * 1000 ;
%
% cell contractility during priming (Pa): stress that cells experience during the extrenal straining period (Pa)
rho_prim_00_24h = rho_prim_00_24h * 1000 ;          
rho_prim_15_24h = rho_prim_15_24h * 1000 ;
rho_prim_30_24h = rho_prim_30_24h * 1000 ;
rho_prim_00_1h  = rho_prim_00_1h  * 1000 ;
rho_prim_15_1h  = rho_prim_15_1h  * 1000 ;
rho_prim_30_1h  = rho_prim_30_1h  * 1000 ; 
%
% durtion of priming (hr)
duration_prim_00_24h = 24 ;
duration_prim_15_24h = 24 ;
duration_prim_30_24h = 24 ;
duration_prim_00_1h  = 1  ;
duration_prim_15_1h  = 1  ;
duration_prim_30_1h  = 1  ;
%
% memory-induced effects: 0% strain for 24h 
tspan = linspace(0,48,10000) ;
[ t_00_24h , y_00_24h ] = ode45(@(t_00_24h,y_00_24h) odefcn( t_00_24h , y_00_24h , duration_prim_00_24h , rho_prim_00_24h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_00_24h = y_00_24h(:,1) ;   
YT_00_24h         = y_00_24h(:,2) ;   
gene1_00_24h      = y_00_24h(:,3) ;    
gene2_00_24h      = y_00_24h(:,4) ;   
%
% memory-induced effects: 15% strain for 24h 
tspan = linspace(0,48,10000) ;
[ t_15_24h , y_15_24h ] = ode45(@(t_15_24h,y_15_24h) odefcn( t_15_24h , y_15_24h , duration_prim_15_24h , rho_prim_15_24h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_15_24h = y_15_24h(:,1) ;   
YT_15_24h         = y_15_24h(:,2) ;   
gene1_15_24h      = y_15_24h(:,3) ;    
gene2_15_24h      = y_15_24h(:,4) ; 
%
% memory-induced effects: 30% strain for 24h 
tspan = linspace(0,48,10000) ;
[ t_30_24h , y_30_24h ] = ode45(@(t_30_24h,y_30_24h) odefcn( t_30_24h , y_30_24h , duration_prim_30_24h , rho_prim_30_24h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_30_24h = y_30_24h(:,1) ;   
YT_30_24h         = y_30_24h(:,2) ;   
gene1_30_24h      = y_30_24h(:,3) ;    
gene2_30_24h      = y_30_24h(:,4) ; 
%
% memory-induced effects: 0% strain for 1h 
tspan = linspace(0,25,10000) ;
[ t_00_1h , y_00_1h ] = ode45(@(t_00_1h,y_00_1h) odefcn( t_00_1h , y_00_1h , duration_prim_00_1h , rho_prim_00_1h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_00_1h = y_00_1h(:,1) ;   
YT_00_1h         = y_00_1h(:,2) ;   
gene1_00_1h      = y_00_1h(:,3) ;    
gene2_00_1h      = y_00_1h(:,4) ;
%
% memory-induced effects: 15% strain for 1h 
tspan = linspace(0,25,10000) ;
[ t_15_1h , y_15_1h ] = ode45(@(t_15_1h,y_15_1h) odefcn( t_15_1h , y_15_1h , duration_prim_15_1h , rho_prim_15_1h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_15_1h = y_15_1h(:,1) ;   
YT_15_1h         = y_15_1h(:,2) ;   
gene1_15_1h      = y_15_1h(:,3) ;    
gene2_15_1h      = y_15_1h(:,4) ;
%
% memory-induced effects: 30% strain for 1h 
tspan = linspace(0,25,10000) ;
[ t_30_1h , y_30_1h ] = ode45(@(t_30_1h,y_30_1h) odefcn( t_30_1h , y_30_1h , duration_prim_30_1h , rho_prim_30_1h , rho_AfterPriming ) , tspan , y0 ) ;
rho_actual_30_1h = y_30_1h(:,1) ;   
YT_30_1h         = y_30_1h(:,2) ;   
gene1_30_1h      = y_30_1h(:,3) ;    
gene2_30_1h      = y_30_1h(:,4) ;
%
% memory-induced effects: baseline 
tspan = linspace(0,24,10000) ;
[ t_base , y_base ] = ode45(@(t_base,y_base) odefcn( t_base , y_base , 24 , rho_AfterPriming , rho_AfterPriming ) , tspan , y0 ) ;
rho_base        = y_base(:,1) ;   
YT_base         = y_base(:,2) ;   
gene1_base      = y_base(:,3) ;    
gene2_base      = y_base(:,4) ;
%
%
% cell activation level 24h after releasing the pre-strain
% see Figure 5C
figure
Y_CellActivation = [rho_actual_00_24h(end) rho_actual_15_24h(end) rho_actual_30_24h(end);rho_actual_00_1h(end) rho_actual_15_1h(end) rho_actual_30_1h(end)] / rho_actual_00_24h(end) ;
bar(Y_CellActivation)
title('Cell activation level')
%
% heat map of cell activation level
CellStress_Range    = linspace(28000,60000,20) ;
CellStress_Duration = linspace(0,25,20)        ;
rho_actual_Matrix            = zeros(length(CellStress_Duration),length(CellStress_Range) ) ;
CellStress_Range_Matrix      = zeros(length(CellStress_Duration),length(CellStress_Range) ) ;
CellStress_Duration_Matrix   = zeros(length(CellStress_Duration),length(CellStress_Range) ) ;
%
for i = 1:length(CellStress_Duration)
    for j = 1:length(CellStress_Range)
        [ tt , yy ] = ode45(@(tt,yy) odefcn( tt , yy , CellStress_Duration(i) , CellStress_Range(j) , rho_AfterPriming ) , linspace(0,CellStress_Duration(i)+24,10000) , y0 ) ;
        rho_actual_Matrix(i,j)          = yy(end,1) / rho_base(end) ;
        CellStress_Range_Matrix(i,j)    = CellStress_Range(j)       ;
        CellStress_Duration_Matrix(i,j) = CellStress_Duration(i)    ;
    end
end
%
% see Figure 4F
figure
contourf( CellStress_Duration_Matrix , CellStress_Range_Matrix/1000 , real(rho_actual_Matrix) , 3 ) ;
xlabel('Cell stress duration (h)')
ylabel('Cell stress, σ_c (kPa)')
%
% relating cell contractility level to tendency of cell to generate contractile force
Constant1 = 5           ;
Constant2 = 0.00000075  ;
%
ActivationLevel_00_24h = ( rho_actual_00_24h(end) - rho_base(end) ) / rho_base(end) ; 
ActivationLevel_15_24h = ( rho_actual_15_24h(end) - rho_base(end) ) / rho_base(end) ; 
ActivationLevel_30_24h = ( rho_actual_30_24h(end) - rho_base(end) ) / rho_base(end) ;
ActivationLevel_00_1h  = ( rho_actual_00_1h(end)  - rho_base(end) ) / rho_base(end) ;
ActivationLevel_15_1h  = ( rho_actual_15_1h(end)  - rho_base(end) ) / rho_base(end) ; 
ActivationLevel_30_1h  = ( rho_actual_30_1h(end)  - rho_base(end) ) / rho_base(end) ;
%
ActivationLevel_00_24h = ( ActivationLevel_00_24h ^ Constant1) * Constant2 ; 
ActivationLevel_15_24h = ( ActivationLevel_15_24h ^ Constant1) * Constant2 ;
ActivationLevel_30_24h = ( ActivationLevel_30_24h ^ Constant1) * Constant2 ; 
ActivationLevel_00_1h  = ( ActivationLevel_00_1h  ^ Constant1) * Constant2 ; 
ActivationLevel_15_1h  = ( ActivationLevel_15_1h  ^ Constant1) * Constant2 ; 
ActivationLevel_30_1h  = ( ActivationLevel_30_1h  ^ Constant1) * Constant2 ; 
% 
% measurement of tissue stress (σ) and cell stress (σ_c) at stationary condition (0% strain) (kPa)
[Stress_00_24h  , CellStress_00_24h ]  = TissueContractileStress( E_M_00         , rho_0 , E_MT , beta*(1-ActivationLevel_00_24h) , alpha ) ; 
[Stress_15_24h  , CellStress_15_24h ]  = TissueContractileStress( E_M_15         , rho_0 , E_MT , beta*(1-ActivationLevel_15_24h) , alpha ) ;
[Stress_30_24h  , CellStress_30_24h ]  = TissueContractileStress( E_M_30         , rho_0 , E_MT , beta*(1-ActivationLevel_30_24h) , alpha ) ; 
[Stress_00_1h   , CellStress_00_1h  ]  = TissueContractileStress( m_beginning_00 , rho_0 , E_MT , beta*(1-ActivationLevel_00_1h)  , alpha ) ; 
[Stress_15_1h   , CellStress_15_1h  ]  = TissueContractileStress( m_beginning_15 , rho_0 , E_MT , beta*(1-ActivationLevel_15_1h)  , alpha ) ;
[Stress_30_1h   , CellStress_30_1h  ]  = TissueContractileStress( m_beginning_30 , rho_0 , E_MT , beta*(1-ActivationLevel_30_1h)  , alpha ) ;
%
% measurement of baseline tissue stress (σ) and cell stress (σ_c) at stationary condition (0% strain) (kPa) 
% the tissue is unconstrained and free to contract without any resistance (σ=0) for the entire period of 24h  
[Stress_Base   , CellStress_Base  ]  = TissueContractileStress( m_beginning_00 , 0.6823*rho_0 , E_MT , beta  , alpha ) ; 
%
% see Figure 1F
figure
Y_TissueStress = 100 * ( [Stress_00_24h Stress_15_24h Stress_30_24h;Stress_00_1h Stress_15_1h Stress_30_1h] - Stress_Base ) / Stress_Base ;
bar(Y_TissueStress)
title('Change in tissue contraction forces (%)')
%
%--------------------------------------------------------------------------
%
% ode function
function dydt = odefun(t,epsilon,E_MT,alpha,beta,epsilon_tot,rho_0, yita, yita1 , yita2 , ksi1 , ksi2 , epsilon2 , transitionPoint)
%
dydt = ((E_MT.*beta-1)./(beta-alpha).*(epsilon_tot-epsilon)+beta.*rho_0/(beta-alpha)-((40*exp(ksi1/yita1*atan(yita1/ksi1))*sqrt(1+(ksi1/yita1)^2)*exp(-ksi1*(epsilon2-epsilon)).*sin(yita1*(epsilon2-epsilon))).*(epsilon<=transitionPoint)+(40*exp(ksi2/yita2*atan(yita2/ksi2))*sqrt(1+(ksi2/yita2)^2)*exp(-ksi2*(epsilon)).*sin(yita2*(epsilon))).*(epsilon>transitionPoint)))./yita;
%
end
%
%--------------------------------------------------------------------------
%
% damage function
function D = DamageFunction( epsilon_M_1 , epsilon_M_2 , D_1 , D_2 , epsilon_M )
%
D = D_1 + ( epsilon_M - epsilon_M_1 ) * ( D_2 - D_1 ) / ( epsilon_M_2 - epsilon_M_1 ) ;
%
end
%
%--------------------------------------------------------------------------
%
% contractile stress in stationary condition
function [S,R] = TissueContractileStress( E_M , rho_0 , E_MT , beta , alpha )
%
%
% tissue force σ
S = E_M * beta * rho_0 / ( E_MT*beta - 1 + E_M*(beta-alpha) ) ;
%
% cell force σ_c
R = ( E_MT + E_M ) * S / E_M ;
%
end
%
%--------------------------------------------------------------------------
%
% signalling netwrok model
% refer to "SI section 2" and equations S11-S14
%
function dydt = odefcn( t , y , duration , rho_prim , rho_current )


dydt=zeros(4,1);

% from under pre-strain to releasing the pre-strain
S=rho_prim.*(t<duration)+rho_current.*(t>=duration) ; 

% cell actomyosin activation level
% denoted as A_c in the maniscript
% see Figure 4D
rho_actual = y(1) ;   
%
% nuclear accumulation of transcriptional and epigenetic factors associated with cell activation
% denoted as T_c in the maniscript
% see Figure 4D
YT         = y(2) ;  
%
% transcription of genes not associated with cell activation
% denoted as g in the maniscript
% see Figure 4D
gene1      = y(3) ;   
%
% transcription of genes associated with cell activation
% denoted as G in the maniscript
% see Figure 4D
gene2      = y(4) ;  
%
% model parameters for the signalling netwrok model
% see Table S2
%
k1 = 0.1 ; k2 = 90 ; k3 = 4 ; k4 = 10 ; k5 = 1 ; 
%
K1 = 1 ; K2 = 11 ; K3 = 100000 ; K4 = 10 ; K5 = 1 ; K6 = 0.5 ; K7 = 160 ; K8 = 45 ;
%
n1 = 5 ; n2 = 1 ; n3 = 5 ; n4 = 138 ; n5 = 5 ; n6 = 8 ; n7 = 3 ; n8 = 1 ;
%
d1 = 0.02 ; d2 = 1 ; d3 = 1 ; d4 = 0.1 ;
%
%
% equations S11-S14 
%
dydt(1)=k1.*((S./K1).^n1+(gene1./K2).^n2)./(1+(S./K1).^n1+(gene1./K2).^n2)+k2.*((S./K3).^n3+(gene2./K4).^n4)./(1+(S./K3).^n3+(gene2./K4).^n4)-d1.*rho_actual;
%
dydt(2)=k3.*rho_actual-d2.*YT;
%
dydt(3)=k4.*(rho_actual./K5).^n5./(1+(rho_actual./K5).^n5+(YT./K6).^n6)-d3.*gene1;
%
dydt(4)=k5.*(YT./K7).^n7./(1+(rho_actual./K8).^n8+(YT./K7).^n7)-d4.*gene2;
%
end

