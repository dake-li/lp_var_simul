function [ABCD_reduced, shock_select] = ABCD_reduce(ABCD)

% Output ABCD representation with reduced dimensions,
% omitting redundant states and shocks

% Starts with ABCD model of the form
        % s_t = A*s_{t-1} + B*e_t
        % y_t = C*s_{t-1} + D*e_t
% with e_t ~ WN(0,I)

% Based on zeros in ABCD matrices, outputs equivalent model
        % x_t = A_reduced*x_{t-1} + B_reduced*epsilon_t
        % y_t = C_reduced*x_{t-1} + D_reduced*epsilon_t
% with fewer states x_t and shocks epsilon_t


A = ABCD.A;
B = ABCD.B;
C = ABCD.C;
D = ABCD.D;

state_select = any(C~=0,1); % States that enter into measurement eqn

% Iteratively expand state vector to include all states that are relevant
% for dynamics of previously selected states
while any(A(state_select,~state_select)~=0,'all') 
    state_select = state_select | any(A(state_select,:)~=0,1);
end
A_reduced = A(state_select,state_select);
B_reduced = B(state_select,:);
C_reduced = C(:,state_select);

% Keep only non-redundant shocks
shock_select = any([B_reduced; D]~=0,1); % Shocks with nonzero loadings
B_reduced = B_reduced(:,shock_select);
D_reduced = D(:,shock_select);

% Output
ABCD_reduced.A = A_reduced;
ABCD_reduced.B = B_reduced;
ABCD_reduced.C = C_reduced;
ABCD_reduced.D = D_reduced;

end