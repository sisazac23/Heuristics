clear
clc



for H=1:20
tic; 
%Lectura de archivos (hojas)
filename = 'Datos.xlsx';
sheet = "I"+num2str(H);
datos = xlsread(filename,sheet);

s = 1000; % Número de soluciones
noise = 50; %Factor de ruido respecto a la distribución uniforme
%Conjunto de soluciones
Sol = [];
for S=1:s
%Lectura de datos
n = datos(1,1); % Número de ítems n
m = datos(1,2); % Número de restricciones m
p = datos(1,3); % Número de funciones objetivo p
x = zeros(1,n); % Candidatos, variable binaria de pertenencia
a = datos(2:m+1,1:n); % Pesos de cada item, para restricción
b = datos(2:m+1,n+1); % Restricción
c = datos(end-(p-1):end,1:n); %Valor de objeto n para función objetivo
b = b';
R = zeros(1,m); %Recursos utilizados por cada restricción
Z = zeros(1,p); %Valores para la función objetivo
UB = zeros(1,p); %Upper Bound - cota superior

 %Vector de soluciones

%Agregar aleatoriedad al problema. Coeficientes de las restricciones
for f=1:m
    for r=1:n
        a(f,r) = a(f,r) + noise*rand(1);
    end
end
%Cálculo de cota superior
for i=1:p
   for j=1:n
       if c(i,j)>0
           UB(1,i) = UB(1,i) + c(i,j);
       end
   end
end

check = zeros(1,n); %Chequea si el elemento n ya se consideró para ingresar o no a la solución
ap = zeros(m,n); %Cambio porcentual entre los elementos disponibles, una vez seleccionado alguno
fact = 0; %Factibilidad de las soluciones

while (sum(check)) < n && (fact==0) 
    Rp = zeros(1,n);
    Rn = zeros(1,n); 
    sel = -1;
    for j=1:n
       if check(1,j)==0
           for i=1:m
               ap(i,j) = a(i,j)/(b(1,i) - R(1,i));
               if (b(1,i) >= 0) && (R(1,i)+a(i,j) > b(1,i))
                  check(1,j) = 1;
                  break
               end
               if (b(1,i) >= 0) &&  (ap(i,j) > Rp(1,j))
                  Rp(1,j) = ap(i,j);
               end
               if (b(1,i) < 0) &&  (ap(i,j) < Rn(1,j))
                  Rn(1,j) = ap(i,j);
               end
           end
           if check(1,j) == 0
               if sel == -1
                   sel = j;
               else
                   if Rn(1,j) > Rn(1,sel)
                       sel = j;
                   else
                       if (Rn(1,j) == Rn(1,sel)) && (Rp(1,j) < Rp(1,sel))
                           sel=j;
                       end
                   end
               end
           end   
       end
    end
    if sel>= 0
        x(1,sel) = 1;
        check(1,sel) = 1;
        fact = 1;
        for i=1:m
            R(1,i) = R(1,i) + a(i,sel);
            if fact == 1 && R(1,i)> b(1,i)
                fact = 0;
            end
        end
    end
end 
%Cálculo de la función objetivo
for k=1:p
    Z(1,k) = sum(c(k,1:end).*x);
end

%Búsqueda de una mejor solución
while sum(check) < n && fact==1
    Rp = ones(1,n).*Inf;
    sel = -1;
    for j=1:n
       if check(1,j) ==0
          for i=1:p
              a(i,j) = (UB(1,i)-Z(1,i))/UB(1,i);
              if Rp(1,j)>ap(i,j)
                  Rp(1,j) = ap(i,j);
              end
          end
          for k=1:m
             if R(1,k)+a(k,j) > b(1,k)
                 check(1,j)=1;
                 break
             end 
          end
          if check(1,j) == 0
              if sel == -1
                  sel = j;
              else
                  if Rp(1,j) > Rp(1,sel)
                      sel = j;
                  end
              end
          end
       end
    end
    if sel >= 0
       x(1,sel) = 1;
       check(sel) = 1;
       fact = 1;
       for i=1:m
          R(1,i) = R(1,i) + a(i,sel);
          if fact == 1 && R(1,i) > b(1,i)
              fact = 0;
          end
       end
       for k=1:p
           Z(1,k) = sum(c(k,1:end).*x);
       end
    end
end
used_r = zeros(1,m);
for k=1:m
   used_r(1,k) = sum(a(k,1:end).*x); 
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
time_fact = [t_fin nfact];
Soluciones(end,1:2) = time_fact;
sheet = "I"+num2str(H);
filename = 'Resultados1';
writematrix(Soluciones,'Resultados1.xls','Sheet', H)
end

