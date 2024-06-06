clear
clc

% Define index of unitary matrix which you want to convert
i_index = 2
j_index = 1

% Define symbolic variables
syms Delta_phi_s_j1 Delta_phi_s_j2 Theta real
syms D1 D2 D3 D4

% Define the basis conversion matrix M_Delta_phi_to_comp
% 00, 01, 10, 11
M_Delta_phi_to_comp = 1/2 * [
    1, -i*exp(i*Delta_phi_s_j2), -i*exp(i*Delta_phi_s_j1), -exp(i*(Delta_phi_s_j1+Delta_phi_s_j2));
    1, i*exp(i*Delta_phi_s_j2), -i*exp(i*Delta_phi_s_j1), exp(i*(Delta_phi_s_j1+Delta_phi_s_j2));
    1, -i*exp(i*Delta_phi_s_j2), i*exp(i*Delta_phi_s_j1), exp(i*(Delta_phi_s_j1+Delta_phi_s_j2));
    1, i*exp(i*Delta_phi_s_j2), i*exp(i*Delta_phi_s_j1), -exp(i*(Delta_phi_s_j1+Delta_phi_s_j2));
];

M_Delta_phi_to_comp = transpose(M_Delta_phi_to_comp);

% Define the U_MS(t) matrix in the Delta_phi basis
U_MS_Delta_phi = [
    D1*exp(i*Theta), 0, 0, 0;
    0, D2*exp(-i*Theta), 0, 0;
    0, 0, D3*exp(-i*Theta), 0;
    0, 0, 0, D4*exp(i*Theta)
];

% Perform the basis transformation
U_MS_comp = simplify(M_Delta_phi_to_comp * U_MS_Delta_phi * ctranspose(M_Delta_phi_to_comp));

% Display the result
disp('The transformed U_MS(t) matrix in the computational basis is:');
disp(U_MS_comp);

% Customize display for D1, D2, D3, D4, Delta_phi_s_j1, Delta_phi_s_j2, Theta
U_MS_comp_latex = latex(4.*U_MS_comp(i_index,j_index));
U_MS_comp_latex = strrep(U_MS_comp_latex, 'D1', 'D_{1}');
U_MS_comp_latex = strrep(U_MS_comp_latex, 'D2', 'D_{2}');
U_MS_comp_latex = strrep(U_MS_comp_latex, 'D3', 'D_{3}');
U_MS_comp_latex = strrep(U_MS_comp_latex, 'D4', 'D_{4}');
%U_MS_comp_latex = strrep(U_MS_comp_latex, 'Delta_phi_s_j1', 'Delta_{phi_{s,j_{1}}}');
%U_MS_comp_latex = strrep(U_MS_comp_latex, 'Delta_phi_s_j2', 'Delta_{phi_{s,j_{2}}}');
U_MS_comp_latex = strrep(U_MS_comp_latex, 'Theta', 'Theta');

% Display the result
disp('The transformed U_MS(t) matrix in the computational basis (LaTeX format) is:');
disp(U_MS_comp_latex);