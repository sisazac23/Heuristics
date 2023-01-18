clear all
clc 

%% Lectura de archivos (hojas)
for H=1:20
filename = 'Datos.xlsx';
sheet = "I"+num2str(H);
datos = xlsread(filename,sheet);
maxiter = 10000;
alpha = 0.15;
s = 1000;
Sol = [];
n = datos(1,1); % Número de ítems n
m = datos(1,2); % Número de restricciones m
p = datos(1,3); % Número de funciones objetivo p
A = datos(2:m+1,1:n); % Pesos de cada item, para restricción
b = datos(2:m+1,n+1); % Restricción
W = datos(end-(p-1):end,1:n); %Valor de objeto n para función objetivo
V = 3; % Número de vecindarios
maxtime = 5; % Tiempo de ejecución, en segundos
t_inicial = tic; % Tiempo de inicio
Z = zeros(1,p); %Valores para la función objetivo
%Criterio para retirar
crit = sum(W);
Cons = zeros(1,m);
final_solutions = []; % Soluciones finales no dominadas
nd_solutions = []; % Soluciones no dominadas de las soluciones de finales
x = ones(1,n);  %Candidatos, variable binaria de pertenencia Todos adentro
totaln = 100;
tnn  = 100;
%% Algoritmo Grasp
while toc(t_inicial) <= maxtime
x = ones(1,n);
Z = zeros(1,p); %Valores para la función objetivo
Cons = zeros(1,m);
while sum(Cons) ~= m
ce = crit.*x;
x_i = find(x);
if alpha > 0
    mincrit = min(ce(x_i));
    maxcrit = max(ce(x_i));
    th = (1-alpha)*(maxcrit - mincrit)+mincrit;
else
    th = max(ce);
end
RCL = (ce<th);
%Mirar elemento para eliminar
if sum(RCL) > 1
    item = find(RCL);
    r = randi([1 length(item)],1,1);
    x(item(r)) = 0;
end
for j=1:m
         v = A(j,1:end).*x;
        if (sum(v)<b(j)) && (b(j)>0) 
            Cons(j) = 1;
        elseif (b(j)<0) && (sum(v)>b(j))
            Cons(j) = 1;
        else
            Cons(j) = 0;
        end
end
maxiter = maxiter-1;
if maxiter==0 || sum(x)<2
    break;
end
end
%Cálculo función objetivo
used_rec = (A*x');
Z = (W*x');
%Cálculo recursos usados
if sum(Cons)==m
   fact = 1;
else
    fact = 0;
end
Sol = [[sum(x)] x used_rec'  Z']; %Solución inicial 
vnd_solutions = [];
final_solutions = [final_solutions; Sol];
vnd_solutions = [vnd_solutions; Sol];
t0 = toc;
j=1;
solutions = 10;
if (H~= 19) && (H~=20)
while solutions >=1
 while j<=V
    if j==1
        [better,x_ls,new_r,new_z] = vnd_add(x,n,W,A,used_rec,b,Z,m);
    elseif j==2
        [better,x_ls,new_r,new_z] = vnd_change(x,n,W,A,used_rec,b,Z,m);
    else
        [better,x_ls,new_r,new_z] = vnd_change2(x,n,W,A,used_rec,b,Z,m);
    end
    if better==true
        x = x_ls;
        used_rec = used_rec + new_r;
        Z = Z + new_z;
        j = 1;
        vnd_solutions = [vnd_solutions; [sum(x_ls)] x_ls used_rec'  Z'];
    else
        j = j+1;
    end
 end
 solutions = solutions -1;
end 
end

[rows,cols] = size(vnd_solutions);
for i=1:rows
    F(i,1:3) = vnd_solutions(i,end-2:end); 
end

ND = paretoDominance(F);
new_S = [];
[row col] = size(ND);
DP = zeros(rows,1);
for r=1:row
    ZF = sum(ND(r,1:end));
    for j=1:rows
    ZS = sum(vnd_solutions(j,end-2:end));
        if ZS == ZF
            DP(j) = 1; 
        end
    end
end
I = find(DP);
new_S = vnd_solutions(I,1:end);
[filas_s cols_s] = size(new_S);
for kk=1:filas_s
final_solutions = [final_solutions; new_S(kk,1:end)];
end
tnn = tnn - 1;
end

[nsol colsol] = size(final_solutions);
x_s = NaN(nsol,n);
for j=1:nsol
x_s(j,1:end) = final_solutions(j,2:n+1);
end
Soluciones = NaN(nsol+1,1+max(final_solutions(1:end,1))+m+p);
Soluciones(1,1) = nsol;
for j=1:nsol
    used_rec1 = final_solutions(j,n+2:end-3);
    FO = final_solutions(j,end-2:end);
    imp = [[sum(x_s(j,1:end))] find(x_s(j,1:end)) used_rec1 FO];
    Soluciones(j+1,1:length(imp)) = imp;
end
sheet = "I"+num2str(H);
writematrix(Soluciones,'vnd_15.xls','Sheet', H);
end

%% Funciónes VND y vecindarios
% 
function [better,x_ls,new_r,new_z] = vnd_add(x,n,W,A,R,b,Z,m)
x0 = x; %Solución inicial
x_ls = x0; %Solución mejorada (si se encuentra una mejor)
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

function [better,x_ls,new_r,new_z] = vnd_change(x0,n,W,A,R,b,Z,m)
x = x0; %Solución inicial
x_ls = x; %Solución mejorada (si se encuentra una mejor)
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

function [better,x_ls,new_r,new_z] = vnd_change2(x0,n,W,A,R,b,z,m)
x = x0; %Solución inicial
x_ls = x; %Solución mejorada (si se encuentra una mejor)
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

%% Pareto dominance
function [ND] = paretoDominance(F)
[i, dim] = size(F);
idxs = [1 : i]';
while i >= 1
    old_size = size(F,1);
    indices = sum( bsxfun( @ge, F(i,:), F ), 2 ) == dim;
    indices(i) = false;
    F(indices,:) = [];
    idxs(indices) = [];
    i = i - 1 - (old_size - size(F,1)) + sum(indices(i:end));
end 
ND = F;
end


