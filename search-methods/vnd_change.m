function [better,x_ls,new_r,new_z] = vnd_change(x0,n,W,A,R,b,Z,m)
x = x0; %Solución inicial
x_ls = x;
z_star = Z; %Función objetivo inicial
new_r = b*0; %Nuevos recursos usados
new_z = Z*0;
% Movimiento posible, se agregan dos elementos y se cambia uno
move = [0 0 0];
% Auxiliar booleano para determinar si se encuentra una mejor solución
better = false;
% Exploración del vecindario de la solución para buscar una mejor
for i = 1:n
    if x(i) == 1
        for j = 1:n
            if x(j) == 0
                for k = 1:n
                    if j ~=k && x(k) == 0
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
                            z_p = Z - W(:,i) + W(:,j) + W(:,k);
                            if prod(z_p>=Z) == 1 && sum(z_p>Z) >= 1
                                better = true;
                                if prod(z_p>=z_star) == 1 && sum(z_p>z_star) >= 1
                                    move(1) = i;
                                    move(2) = j;
                                    move(3) = k;
                                    z_star = z_p;
                                end
                            end
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
    k = move(3);
    x(i) = 0; %Saco un elemento de la bolsa
    x(j) = 1; %Agrego un elemento a la bolsa
    x(k) = 1; %Agrego un elemento a la bolsa
    x_ls = x;
    new_r = A(:,j) + A(:,k) - A(:,i);
    new_z = W(:,j) + W(:,k) - W(:,i);
end
end