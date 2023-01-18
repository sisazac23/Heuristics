function [better,x_ls,new_r,new_z] = vnd_add(x,n,W,A,R,b,Z,m)
x0 = x; %Solución inicial
x_ls = x0;
z_ini = Z; % Función objetivo inicial
new_r = b*0; %Recursos usados por el nuevo item
new_z = Z*0; % Peso en función objetivo
% Movimiento posible 
move = 0;
% Auxiiar booleano, verdadero si se encuentra una solución
better = false;
% Exploración del vecindario de la solución para buscar una mejor
for i = 1:n
    if x(i) == 0
        check = zeros(m,1);
        for j=1:m
            if ((R(j) + A(j,i) <= b(j)) && (b(j)>0))
                check(j,1) = 1;
            elseif ((R(j) + A(j,i) >= b(j)) && (b(j)<0))
                check(j,1) = 1;
            else
                check(j,1) = 0;
            end
        end
        if sum(check) == m
            z_p = Z + W(:,j);
            if prod(z_p>=Z) == 1 && sum(z_p>Z) >= 1
                better = true;
                if prod(z_p>=z_ini) == 1 && sum(z_p>z_ini) >= 1
                    move = i;
                    z_ini = z_p;
                end
            end
        end
    end
end
% Si se encuentra una mejor solución (que sea dominante)
if better == true
    i = move;
    x0(i) = 1;
    x_ls = x0;
    new_r = A(:,i);
    new_z = W(:,i);
end
end

