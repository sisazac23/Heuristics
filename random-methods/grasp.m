clear
clc


for H=1:20
tic;
%Lectura de archivos (hojas)
filename = 'Datos.xlsx';
sheet = "I"+num2str(H);
datos = xlsread(filename,sheet);
maxiter = 10000;
alpha = 0.25;
s = 1000;
Sol = [];
n = datos(1,1); % Número de ítems n
m = datos(1,2); % Número de restricciones m
p = datos(1,3); % Número de funciones objetivo p
a = datos(2:m+1,1:n); % Pesos de cada item, para restricción
b = datos(2:m+1,n+1); % Restricción
c = datos(end-(p-1):end,1:n); %Valor de objeto n para función objetivo
for S=1:s
R = zeros(1,m); %Recursos utilizados por cada restricción
Z = zeros(1,p); %Valores para la función objetivo
%Criterio para retirar
crit = sum(c);
Cons = zeros(1,m);
x = ones(1,n);  %Candidatos, variable binaria de pertenencia Todos adentro
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
         v = a(j,1:end).*x;
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

for i=1:p
    Z(i) = sum(c(i,1:end).*x);
end
used_r = zeros(1,m);
for k=1:m
   used_r(1,k) = sum(a(k,1:end).*x); 
end
if sum(Cons)==m
   fact = 1;
else
    fact = 0;
end
    Sol = [Sol; [sum(x)] x used_r  Z fact];
end
%%Soluciones factibles o infactibles
it = 1;
for k=1:s
    if sum(Sol(1:end,end)) == 0
       F = Sol(1:end,end-3:end-1); 
       break;
    end
    if Sol(k,end) == 1
        F(it,1:3) = Sol(k,end-3:end-1);
        it = it+1;
    end
end
nfact = sum(Sol(1:end,end)); %Número de soluciones factibles
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


new_S = [];
[row col] = size(F);
DP = zeros(s,1);
for r=1:row
    ZF = sum(F(r,1:end));
    for j=1:s
    ZS = sum(Sol(j,end-3:end-1));
        if ZS == ZF
            DP(j) = 1; 
        end
    end
end

I = find(DP);
new_S = Sol(I,1:end-1);
x_s = NaN(row,n);
for j=1:row
x_s(j,1:end) = new_S(j,2:n+1);
end
%Escribir las respuestas
%Número de soluciones
nsol = row;
%Formato de solución
Soluciones = NaN(nsol+1+1,1+max(new_S(1:end,1))+m+p);
Soluciones(1,1) = nsol;
for j=1:nsol
    rec_usados = new_S(j,n+2:end-3);
    FO = new_S(j,end-2:end);
    imp = [[sum(x_s(j,1:end))] find(x_s(j,1:end)) rec_usados FO];
    Soluciones(j+1,1:length(imp)) = imp;
end
toc;
t_fin=toc;
time_fact = [t_fin nfact alpha];
Soluciones(end,1:3) = time_fact;
sheet = "I"+num2str(H);
filename = 'grasp_25';
writematrix(Soluciones,'grasp_25.xls','Sheet', H)
end
