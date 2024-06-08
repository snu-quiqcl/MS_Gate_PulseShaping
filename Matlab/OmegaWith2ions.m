clear
clc

% Number of m1, m2 list

N = 6;
LIST_NUM = N;
DISPLAY = "delta";

% Define m1_list and m2_list
m1_list = 1:1:N;
m2_list = 1:1:N;
m1 = 6;
m2 = 1;

% Define constants (replace these with actual values)
omega_1 = 2 .* pi .* 2.05 .* 1E6;
omega_2 = 2 .* pi .* 2.132 .* 1E6;
M = 170.936323 .* (1E-3)./ 6.02E23;
delta_k = 28339146.473469555;
hbar = 6.626E-34./(2.*pi);

% Calculate Omega values using array operations
%delta = (M2 ./ (M1 + M2)) * omega_1 + (M1 ./ (M1 + M2)) * omega_2;
%tau = 2 * (M1 + M2) / (omega_1 - omega_2);
%Omega = sqrt((pi * omega_1.*omega_2 * M ./ (abs(omega_2 - omega_1) * abs(delta .* tau - sin(delta .* tau)) * delta_k^2 * hbar))) .* delta;

delta = zeros(LIST_NUM,LIST_NUM);
tau = zeros(LIST_NUM,LIST_NUM);
Omega = zeros(LIST_NUM,LIST_NUM);

for i = 1:LIST_NUM
    for j = 1:LIST_NUM
        if (m1_list(i) + m2_list(j)) == 0
            delta(i,j) = NaN;
            tau(i,j) = NaN;
            Omega(i,j) = NaN;
        else
            delta(i,j) = ( ...
                ( ...
                    m2_list(j) ./ ( m1_list(i) + m2_list(j) ) ...
                ) .* omega_1 ...
                + ...
                ( ...
                    m1_list(i) ./ ( m1_list(i) + m2_list(j) ) ...
                ) .* omega_2 ...
            );
            
            tau(i,j) = 2 .* pi .* abs((m1_list(i) + m2_list(j)) ./ (omega_1 - omega_2));
            
            A = 0;
            A = 1./( omega_1 .* ( omega_1 - delta(i,j) ) ) - 1./( omega_2 .* ( omega_2 - delta(i,j) ) );
            
            Omega(i,j) = sqrt(2.*M.*pi ./ (hbar .*tau(i,j).* abs(A) ) ) ./ delta_k;
        end
    end
end

% Create 3D bar plot
if strcmp(DISPLAY,"Omega") == 1
    figure;
    disp(min(min(Omega./(pi * 2 * 1E3))));
    b = bar3(Omega./(pi * 2 * 1E3));
    set(gca,'XTickLabel',m1_list) 
    set(gca,'YTickLabel',m2_list) 
    
    % Set axis labels and title
    xlabel('m1');
    ylabel('m2');
    zlabel('\Omega [kHz]');
    title('3D Bar Plot of \Omega');
elseif strcmp(DISPLAY,"tau") == 1
    figure;
    disp(max(max(tau.* 1E6)));
    b = bar3(tau.* 1E6);
    set(gca,'XTickLabel',m1_list) 
    set(gca,'YTickLabel',m2_list) 
    
    % Set axis labels and title
    xlabel('m1');
    ylabel('m2');
    zlabel('\tau [\mu s]');
    title('3D Bar Plot of \tau');
elseif strcmp(DISPLAY,"delta") == 1
    figure;
    disp(min(min(delta./(pi * 2 * 1E6))));
    b = bar3(delta./(pi * 2 * 1E6));
    set(gca,'XTickLabel',m1_list) 
    set(gca,'YTickLabel',m2_list) 
    
    % Set axis labels and title
    xlabel('m1');
    ylabel('m2');
    zlabel('\delta [MHz]');
    zlim([(omega_1 .* 0.999 ./(2 * pi * 1E6)) ( omega_2 .* 1.001  ./(2 * pi * 1E6))])
    title('3D Bar Plot of \delta');
end

% Set color based on Omega values
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% Adjust the colorbar
colorbar;