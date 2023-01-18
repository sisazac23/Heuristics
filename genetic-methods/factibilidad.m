function fact=factibilidad(used_rec,b,m)
% Revisa la factibilidad de una solución
% Retorna un valor booleano, 1 si es factible, 0 si no
Cons = zeros(1,m);
b = b';
used_rec = used_rec';
for i=1:m
    if b(1,i)>0
        if used_rec(1,i)<=b(1,i)
            Cons(1,i)= 1;
        end
    end
    if b(1,i)< 0
       if used_rec(1,i)>= b(1,i)
           Cons(1,i) = 1;
       end
    end
if sum(Cons)==m
   fact = 1;
else
    fact = 0;
end
end