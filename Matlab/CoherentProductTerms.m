clear
clc

%% Calculate coherent state inter product
syms a1 a2

A1 =  a1 + a2;
A2 =  a1 - a2;
A3 = -a1 + a2;
A4 = -a1 - a2;

A_list = [A1 A2 A3 A4];

for i = 1:4
    for j = 1:4
        A = A_list(i);
        B = A_list(j);
        
        A_B = simplify(- conj(A) * A ./ 2 - conj(B) * B ./ 2 + conj(A) * B);
        disp(" ")
        switch i
            case 1
                disp("\braket{\alpha_{1}+\alpha_{2}|")
            case 2
                disp("\braket{\alpha_{1}-\alpha_{2}|")
            case 3
                disp("\braket{-\alpha_{1}+\alpha_{2}|")
            case 4
                disp("\braket{-\alpha_{1}-\alpha_{2}|")
        end

        switch j
            case 1
                disp("\alpha_{1}+\alpha_{2}}")
            case 2
                disp("\alpha_{1}-\alpha_{2}}")
            case 3
                disp("-\alpha_{1}+\alpha_{2}}")
            case 4
                disp("-\alpha_{1}-\alpha_{2}}")
        end
        A_B = strrep(latex(A_B),"a_{1}","\alpha_{1}");
        A_B = strrep(A_B,"a_{2}","\alpha_{2}");
        
        disp(" = ")
        disp("\exp(" + string(A_B) + ")")
        disp("\\")
    end
end