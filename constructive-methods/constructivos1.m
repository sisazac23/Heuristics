clear
clc


for H=1:20
tic; 
filename = 'Datos.xlsx';
sheet = "I"+num2str(H);
datos = xlsread(filename,sheet);

n = datos(1,1); % Número de ítems n
m = datos(1,2); % Número de restricciones m
p = datos(1,3); % Número de funciones objetivo p
x_i = ones(1,n); % Candidatos, variable binaria de pertenencia
aik = datos(2:m+1,1:n); % Pesos de cada item, para restricción
bk = datos(2:m+1,n+1); % Restricción
wih = datos(end-(p-1):end,1:n); %Valor de objeto n para función objetivo

R = zeros(1,length(bk)); %Vector binario con 0 si no cumple restricción, 1 si sí
KK = sum(wih); %Suma de los pesos de Wih para cada item i
[S,I] =sort(KK(1,1:end),'descend'); %Criterio de selección: según peso en función objetivo
i = n; %Contador para recorrer los objetos
while sum(R)~= length(bk)
    %x_i(I(i)) = 0;
    for j=1:length(bk)
            v = aik(j,1:end).*x_i;
        if (sum(v)<bk(j)) && (bk(j)>0) 
            R(j) = 1;
            %x_i(I(i)) = 1;
            
        elseif (bk(j)<0) && (sum(v)>bk(j))
            R(j) = 1;
            %x_i(I(i)) = 1;
            
        else
            R(j) = 0;
            x_i(I(i))=0;
            
        end
    end
    i = i - 1;
    if i==1
       break
    end
end
Z = zeros(1,p);
for i=1:p
    Z(i) = sum(wih(i,1:end).*x_i);
end
toc;
t_fin=toc;
%%Ahora escribamos estos resultados en un archivo de Excel
A = [1];
sheet="Hoja"+num2str(H);
xlswrite("Resultados1",A,sheet,'A1'); %Escribimos número de soluciones

B = [[sum(x_i)]];
xlswrite("Resultados1",B,sheet,'A2'); %Escribimos número de elementos usados

C = [find(x_i==1)];
xlswrite("Resultados1",C,sheet,"A3"); %Escribimos índice de elementos usados

Recursos = [];
for r=1:length(bk)
   Recursos(r) = sum(aik(r,1:end).*x_i); 
end
xlswrite("Resultados1",Recursos,sheet,'A4') %Escribimos la cantidad de recursos usados por restricción

FO = [Z];
xlswrite("Resultados1",FO,sheet,"A5"); %Escribimos los valores de la función objetivo

xlswrite("Resultados1",t_fin,sheet,"A6"); %Escribimos el tiempo que tardó en resolver el problema
end

