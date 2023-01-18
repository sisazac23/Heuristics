function [padre1,padre2] = torneo(X,Rank)
% Proceso de selección de las soluciones padres
% Una vez que las soluciones estén clasificadas se realizarán permutaciones
% para escoger parejas de ellas. Donde en el proceso se busca escoger a
% aquellas cuya suma de ranking sea lo menor posible (lo más cercano a 1)
p = size(X,1);
rp = randperm(p,2);
sp = randperm(p,2);
r1 = sum(Rank(rp'));
r2 = sum(Rank(sp'));
if r1 < r2
    padre1 = X(rp(1),:);
    padre2 = X(rp(2),:);
else 
    padre1 = X(sp(1),:);
    padre2 = X(sp(2),:);
end
end