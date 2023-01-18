function [better,x_ls,new_r,new_z] = vnd_change2(x0,n,W,A,R,b,z,m)
x = x0; %Solución inicial
x_ls = x;
z_star = z; % Función objetivo inicial
new_r = b*0; % Nuevos recursos usados
new_z = z*0; % Nueva función objetivo
move = [0 0]; % Movimientos a realizar
better = false; % Variable auxiliar booleana que determina si se halló o no una mejor solución
for i = 1:n
    if x(i) == 1
        for j = 1:n
            if x(j) == 0
                check = zeros(m,1);
                for u=1:m
                    if ((R(u) + A(u,j) <= b(u)) && (b(u)>0))
                        check(u,1) = 1;
                    elseif ((R(u) + A(u,j) >= b(u)) && (b(u)<0))
                        check(u,1) = 1;
                    else
                        check(u,1) = 0;
                   end
                end
                if sum(check) == m
                    z_p = z + W(:,j) - W(:,i);
                    if prod(z_p>=z) == 1 && sum(z_p>z) >= 1
                        better = true;
                        if prod(z_p>=z_star) == 1 && sum(z_p>z_star) >= 1
                            move(1) = i;
                            move(2) = j;
                            z_star = z_p;
                        end
                    end
                end
            end
        end
    end
end
if better == true
    i = move(1);
    j = move(2);
    x(i) = 0; % Elemento que saco de la bolsa
    x(j) = 1; % Elemento que agrego a la bolsa
    x_ls = x;
    new_r = A(:,j) - A(:,i);
    new_z = W(:,j) - W(:,i);
end
end