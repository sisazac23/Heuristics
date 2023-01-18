function x = mutacion(x,p,n)
% Operador genético de mutación
% Cambia el valor de un item si se cumple el valor condicional relaciona con
% el parámetro p_mutacion, 0 si antes el item estaba seleccionado y un 1 si antes ese item no
% estaba.
for i = 1:n
    r = rand;
    if r <= p
        x(i) = ~x(i);
    end
end
end