%Symmetric

syms mu_c cbar V C_z_alphadot C_m_alphadot K_y
C1_sym = [  -2*mu_c*cbar/V^2        0                               0               0;...
            0                       (C_z_alphadot-2*mu_c)*cbar/V    0               0;...
            0                       0                               -cbar/V         0;...
            0                       C_m_alphadot*cbar/V             0               -2*mu_c*K_y^2*cbar^2/V^2];
        

syms C_X_u C_X_alpha C_Z_0 C_X_q C_Z_u C_Z_q C_Z_alpha C_X_0 C_m_u C_m_alpha C_m_q

C2_sym = [C_X_u     C_X_alpha   C_Z_0       C_X_q*cbar/V;...
            C_Z_u   C_z_alphadot    -C_X_0  (C_Z_q+2*mu_c)*cbar/V; ...
            0       0               0       cbar/V;...
            C_m_u   C_m_alpha   0   C_m_q*cbar/V];
        
A= inv(C1_sym)*C2_sym;
            