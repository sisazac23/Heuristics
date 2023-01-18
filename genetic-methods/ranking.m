function paretoF = ranking(Z)
% Ranking de soluciones en frentes de Pareto
% Dominancia de Pareto para hallar soluciones no dominadas
% Luego estas reciben el primer rango, pues pertenecen al primer frente de
% Pareto. Luego a los valores que ya se clasificaron se les dan valores muy pequeños para que no sean tenidos
% para futuros frentes, la clasificación termina cuando todos los items
% tengan un valor de clasificación.
[k,p] = size(Z);
paretoF = zeros(k,1);
PD = Z;
contador = 1;
rango = 0;
ND = paretoDominance(PD);
[fila,~] = size(ND);
DP = zeros(k,1);
ZS = sum(Z(:,1:end-1),2);
for r=1:fila
    ZF = sum(ND(r,1:end-1));
    for j=1:k
        if ZS(j) == ZF
            DP(j) = 1; 
        end
    end
end
I = find(DP);
paretoF(I) = contador;
PD(I,:) = (-100000)*ones(length(I),p);
contador = contador + 1;
rango = rango + length(I);
while rango ~= k
    ND = paretoDominance(PD);
    [fila,~] = size(ND);
    DP = zeros(k,1);
    ZS = sum(Z(:,1:end-1),2);
    for r=1:fila
        ZF = sum(ND(r,1:end-1));
        for j=1:k
            if ZS(j) == ZF
                DP(j) = 1;
            end
        end
    end
    I = find(DP);
    paretoF(I) = contador;
    PD(I,:) = (-100000)*ones(length(I),p);
    contador = contador + 1;
    rango = rango + length(I);
end
end