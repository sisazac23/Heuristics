function hijo = cruce(padre1,padre2,p_cruce,n,m,W,A,b,p,alpha,maxiter)
% El cruce se hace en un solo punto (single-point crossver)
% Las soluciones generadas no necesariamente son factibles
% Adicionalmente se usa la distancia de Hamming para conocer la cercanía
% entre dos arreglos binarios (los dos padres), si la distancia no cumple el criterio,
% el hijo será producto de un métood aleatorizado (GRASP)
hijo = zeros(1,n);
distanciag = pdist([double(padre1); double(padre2)],'hamming');
if distanciag > p_cruce
    % Punto aleatorio de cruce
    pc = randi(n-1);
    % Cruce de un punto
    hijo(1:pc) = padre1(1:pc);
    hijo(pc+1:end) = padre2(pc+1:end);
else
    % Nuevo hijo aleatorio 
     [hijo,~,~,~] = grasp(W,A,b,n,p,m,maxiter,alpha);
end

end