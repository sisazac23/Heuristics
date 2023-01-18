clear
clc

for H=1:20
filename = 'Datos.xlsx';
sheet = "I"+num2str(H);
datos = xlsread(filename,sheet);
maxiter = 10000;
alpha = 0.20;
s = 1000;
Sol = [];
n = datos(1,1); % Número de ítems n
m = datos(1,2); % Número de restricciones m
p = datos(1,3); % Número de funciones objetivo p
A = datos(2:m+1,1:n); % Pesos de cada item, para restricción
b = datos(2:m+1,n+1); % Restricción
W = datos(end-(p-1):end,1:n); %Valor de objeto n para función objetivo
V = 3; % Número de vecindarios
maxtime = 3; % Tiempo de ejecución, en segundos
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
% Perturbación
solutions = 10;
j=1;
if (H~= 19) && (H~=20)
%j=1;
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
 solutions = solutions - 1;
end 
end
[nper colper] = size(vnd_solutions);

idx = randperm(n,4);
x_per = vnd_solutions(1,2:n+1);
for i = idx
   x_per(i) = ~x_per(i); 
end
%Cálculo función objetivo
used_rec2 = (A*x_per');
Z2 = (W*x_per');
solutions=10;
vnd_solutions2 = vnd_solutions;
k=1;
if (H~= 19) && (H~=20)
while solutions >=1
 %k=1;
 while k<=V
    if k==1
        [better,x_ls,new_r,new_z] = vnd_add(x_per,n,W,A,used_rec2,b,Z2,m);
    elseif k==2
        [better,x_ls,new_r,new_z] = vnd_change(x_per,n,W,A,used_rec2,b,Z2,m);
    else
        [better,x_ls,new_r,new_z] = vnd_change2(x_per,n,W,A,used_rec2,b,Z2,m);
    end
    if better==true
        x2 = x_ls;
        used_rec2 = used_rec2 + new_r;
        Z2 = Z2 + new_z;
        k = 1;
        vnd_solutions2 = [vnd_solutions2; [sum(x2)] x2 used_rec2'  Z2'];
    else
        k = k+1;
    end
 end
 solutions = solutions - 1;
end 
end
[rows,cols] = size(vnd_solutions2);
for i=1:rows
    F(i,1:3) = vnd_solutions2(i,end-2:end); 
end

ND = paretoDominance(F);
new_S = [];
[row col] = size(ND);
DP = zeros(rows,1);
for r=1:row
    ZF = sum(ND(r,1:end));
    for j=1:rows
    ZS = sum(vnd_solutions2(j,end-2:end));
        if ZS == ZF
            DP(j) = 1; 
        end
    end
end
I = find(DP);
new_S = vnd_solutions2(I,1:end);
[filas_s cols_s] = size(new_S);
for kk=1:filas_s
final_solutions = [final_solutions; new_S(kk,1:end)];
end
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
writematrix(Soluciones,'ils_20.xls','Sheet', H);
end
