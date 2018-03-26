%% Symmetrical - Generate A and B matrix from C1 - C2 - C3

C1_sym = [  -2*muc*c/V0^2        0                               0               0;...
            0                       (CZadot-2*muc)*c/V0    0               0;...
            0                       0                               -c/V0         0;...
            0                       Cmadot*c/V0             0               -2*muc*KY2*c^2/V0^2];
        
        
C2_sym = [CXu/V0     CXa   CZ0       CXq*c/V0;...
            CZu/V0   CZa    -CX0  (CZq+2*muc)*c/V0; ...
            0       0               0       c/V0;...
            Cmu/V0   Cma   0   Cmq*c/V0];
        
C3_sym = [CXde ; CZde ; 0 ; Cmde ];

% Calculate A and B matrices

A = -inv(C1_sym) * C2_sym;
B = -inv(C1_sym) * C3_sym;

 C = eye(4);
 D = zeros(4,1);