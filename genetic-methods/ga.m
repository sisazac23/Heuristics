clear all
clc

%% Lectura de datos
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
%% Valores y parámetros
V = 3; % Número de vecindarios
maxtime = 300; % Tiempo de ejecución, en segundos
maxiter = 1000;
t_inicial = tic; % Tiempo de inicio
generacion = 100;
poblacion = 200;
num_hijos = 50;
p_cruce = 0.15; %Porcentaje de cruce entre dos padres, si la cercanía genética es menor a este valor no se reproducirán
epsilon = 0.5; 
p_mutacion = 0.05; %Probabilidad de mutación
X = zeros(poblacion,n);
Zf = zeros(poblacion,p+1);
nsol = 0; % Número de soluciones
Xgen = zeros(poblacion,n);
Zgen = zeros(poblacion,p+1);
%Generar población inicial - método grasp
t_inicial = tic; % Tiempo de inicio
while toc(t_inicial) <= maxtime
    for i =1:poblacion
        [x,used_rec,Z,fact] = grasp(W,A,b,n,p,m,maxiter,alpha);
        Xgen(i,:) = x;
        Zgen(i,:) = [Z' fact];
        nsol = nsol + 1;
    end
    if H>=19
        X = Xgen;
        Zf = Zgen;
    end
    if H <19
        Rank = ranking(Zgen);
        for gen = 1:generacion
            if std(Rank)<epsilon || toc(t_inicial) > maxtime
                break;
            end
            % Hijos
            descendencia = zeros(num_hijos,n);
            Zson = zeros(num_hijos,p+1);
            for i = 1:num_hijos
                % Torneo para seleccionar padres
                [padre1,padre2] = torneo(Xgen,Rank);
                hijo = cruce(padre1,padre2,p_cruce,n,m,W,A,b,p,alpha,maxiter);
                fact_hijo = factibilidad(A*hijo',b,m);
                if fact_hijo == 1
                    % Mejora de la solución
                    neigh = 1;
                    while neigh<=V
                        if neigh==1
                            [better,x_hijo,~,~] = vnd_add(hijo,n,W,A,A*hijo',b,W*hijo',m);
                        elseif neigh==2
                            [better,x_hijo,~,~] = vnd_change(hijo,n,W,A,A*hijo',b,W*hijo',m);
                        else
                            [better,x_hijo,~,~] = vnd_change2(hijo,n,W,A,A*hijo',b,W*hijo',m);
                        end
                        if better==true
                            fact_hijo = factibilidad(A*x_hijo',b,m);
                            if fact_hijo==1
                                hijo = x_hijo;
                                break;
                            end
                        else
                            neigh = neigh+1;
                        end
                    end
                end
                % Posible mutación genética
                r = randn;
                if r <= p_mutacion
                    hijo = mutacion(hijo,p_mutacion,n);
                    fact_hijo = factibilidad(A*hijo',b,m);
                end
                descendencia(i,:) = hijo;
                Zson(i,:) = [(W*hijo')' fact_hijo];
                nsol = nsol + 1;
            end
            Xgen = [Xgen; descendencia];
            Zgen = [Zgen; Zson];
            infact = find(Zgen(:,end)~=1);
            [~,sizeZ] = size(Zgen);
            Zgen(infact,:) = (-100000)*ones(length(infact),sizeZ);
            Rank = ranking(Zgen);
            [~, newPoblacion] = sort(Rank);
            newPoblacion = newPoblacion(1:poblacion);
            Xgen = Xgen(newPoblacion,:);
            Zgen = Zgen(newPoblacion,:);
            Rank = Rank(newPoblacion,:);
            if toc(t_inicial) > maxtime
                break;
            end
        end
        
        %% Mejora de las soluciones
        for i=1:poblacion
            x = Xgen(i,:);
            fact = Xgen(i,p+1);
            if fact == 1
                neigh = 1;
                while neigh<=V
                    if neigh==1
                        [better,x_ls,~,~] = vnd_add(x,n,W,A,A*x',b,W*x',m);
                    elseif neigh==2
                        [better,x_ls,~,~] = vnd_change(x,n,W,A,A*x',b,W*x',m);
                    else
                        [better,x_ls,~,~] = vnd_change2(x,n,W,A,A*x',b,W*x',m);
                    end
                    if better==true
                        fact_ls = factibilidad(A*x_ls',b,m);
                        if fact_ls==1
                            x = x_ls;
                            break;
                        end
                    else
                        neigh = neigh+1;
                    end
                end
            end
            fact_ls = factibilidad(A*x',b,m);
            Xgen = [Xgen; x];
            Zgen = [Zgen; (W*x')' fact_ls];
        end
        X = [X; Xgen];
        Zf = [Zf; Zgen];
    end
end
%% Dominancia de Pareto
[X,ix,~] = unique(X,'rows');
Zf = Zf(ix,:);
[sizeZ,~] = size(Zf);
new_S = [];
ND = paretoDominance(Zf);
[row,~] = size(ND);
DP = zeros(sizeZ,1);
ZS = sum(Zf(:,1:end-1),2);
for r=1:row
    ZF = sum(ND(r,1:end-1));
    for j=1:sizeZ
        if ZS(j) == ZF
            DP(j) = 1; 
        end
    end
end
%% Escritura de resultados
I = find(DP);
new_S = X(I,1:end);
Soluciones = NaN(row+1,1+max(sum(new_S,2))+m+p);
Soluciones(1,1) = row;
for j=1:row
    imp = [[sum(new_S(j,1:end))] find(new_S(j,1:end)) (A*new_S(j,1:end)')'  (W*new_S(j,1:end)')'];
    Soluciones(j+1,1:length(imp)) = imp;
end
sheet = "I"+num2str(H);
writematrix(Soluciones,'ga05.xls','Sheet', H)
end

