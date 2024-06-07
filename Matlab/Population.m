clear
clc

%% Calculate polpulation with braket terms
% P(i_index, j_index)
i_index = 0
j_index = 1

syms A1 A2 A3 A4
% A1 =  a1 + a2;
% A2 =  a1 - a2;
% A3 = -a1 + a2;
% A4 = -a1 - a2;
syms Theta Delta_phi_s_j1 Delta_phi_s_j2 real

% P00
if i_index == 0 && j_index == 0
    x = exp(-1i*Theta) .* A2 + exp(-1i*Theta) .* A3 + exp(1i*Theta) .* A1 + exp(1i*Theta) .* A4;
end

%P10
if i_index == 1 && j_index == 0
    x = ( - exp(1i*Delta_phi_s_j1) .* exp(1i*Theta) .* A1 ...
        - exp(1i*Delta_phi_s_j1) .* exp(-1i*Theta) .* A2 ...
        + exp(1i*Delta_phi_s_j1) .* exp(-1i*Theta) .* A3 ...
        + exp(1i*Delta_phi_s_j1) .* exp(1i*Theta) .* A4 );
end

%P01
if i_index == 0 && j_index == 1
    x = ( - exp(1i*Delta_phi_s_j1) .* exp(1i*Theta) .* A1 ...
        + exp(1i*Delta_phi_s_j1) .* exp(-1i*Theta) .* A2 ...
        - exp(1i*Delta_phi_s_j1) .* exp(-1i*Theta) .* A3 ...
        + exp(1i*Delta_phi_s_j1) .* exp(1i*Theta) .* A4 );
end

% P11
if i_index == 1 && j_index == 1
    x = ( -exp(1i*Delta_phi_s_j1) .* exp(1i*Delta_phi_s_j2) .* exp(1i*Theta) .* A1 ...
        + exp(1i*Delta_phi_s_j1) .* exp(1i*Delta_phi_s_j2) .* exp(-1i*Theta) .* A2 ...
        + exp(1i*Delta_phi_s_j1) .* exp(1i*Delta_phi_s_j2) .* exp(-1i*Theta) .* A3 ... 
        - exp(1i*Delta_phi_s_j1) .* exp(1i*Delta_phi_s_j2) .* exp(1i*Theta) .* A4 );
end

S = expand(conj(x) * x);
%simplify(expand(conj(x) * x),'Steps',100,'All',true)
S = latex(S);
S = strrep(S,"\overline{A_{1}}","\bra{\alpha_{1}+\alpha_{2}}");
S = strrep(S,"\overline{A_{2}}","\bra{\alpha_{1}-\alpha_{2}}");
S = strrep(S,"\overline{A_{3}}","\bra{-\alpha_{1}+\alpha_{2}}");
S = strrep(S,"\overline{A_{4}}","\bra{-\alpha_{1}-\alpha_{2}}");
S = strrep(S,"A_{1}","\ket{\alpha_{1}+\alpha_{2}}");
S = strrep(S,"A_{2}","\ket{\alpha_{1}-\alpha_{2}}");
S = strrep(S,"A_{3}","\ket{-\alpha_{1}+\alpha_{2}}");
S = strrep(S,"A_{4}","\ket{-\alpha_{1}-\alpha_{2}}");
% Swap conj(A) * A to A * conj(A) even with other patterns in between
fun = @(s)sprintf('.{%d}',find(cumsum((s=='}' )-(s=='{'))>0,1,'first'));
S = regexprep(S, '\\ket\{((??@fun($'')))(.+?(?=\\bra))\\bra\{((??@fun($'')))', '$2\\braket{{$3\|{$1}');
S = regexprep(S, '\\bra\{((??@fun($'')))(.+?(?=\\ket))\\ket\{((??@fun($'')))', '$2\\braket{{$1\|{$3}');

% with braket term
disp("Result in braket Terms")
disp(S)

%% Convert braket terms to explicit terms

syms a1 a2

A_list = [a1 + a2, a1 - a2, -a1 + a2 ,-a1 - a2];
braket_list = strings(4:4);
braket_map_list = strings(4:4);

for i = 1:4
    for j = 1:4
        A = A_list(i);
        B = A_list(j);
        
        A_B = simplify(- conj(A) * A ./ 2 - conj(B) * B ./ 2 + conj(A) * B);
        switch i
            case 1
                braket_list(i,j) = "\braket{{\alpha_{1}+\alpha_{2}}|";
            case 2
                braket_list(i,j) = "\braket{{\alpha_{1}-\alpha_{2}}|";
            case 3
                braket_list(i,j) = "\braket{{-\alpha_{1}+\alpha_{2}}|";
            case 4
                braket_list(i,j) = "\braket{{-\alpha_{1}-\alpha_{2}}|";
        end

        switch j
            case 1
                braket_list(i,j) = braket_list(i,j) + "{\alpha_{1}+\alpha_{2}}}";
            case 2
                braket_list(i,j) = braket_list(i,j) + "{\alpha_{1}-\alpha_{2}}}";
            case 3
                braket_list(i,j) = braket_list(i,j) + "{-\alpha_{1}+\alpha_{2}}}";
            case 4
                braket_list(i,j) = braket_list(i,j) + "{-\alpha_{1}-\alpha_{2}}}";
        end
        A_B = strrep(latex(A_B),"a_{1}","\alpha_{1}");
        A_B = strrep(A_B,"a_{2}","\alpha_{2}");
        
        braket_map_list(i,j) = ("\exp(" + string(A_B) + ")\\");
    end
end

for i = 1:4
    for j = 1:4
        S = strrep(S,braket_list(i,j),braket_map_list(i,j));
    end
end

disp("Result in explicit Terms")
disp(S)