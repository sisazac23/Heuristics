function [x,used_rec,Z,fact] = grasp(W,A,b,n,p,m,maxiter,alpha)
x1 = ones(1,n);
crit = sum(W);
Cons = zeros(1,m);
Z = zeros(1,p); %Valores para la función objetivo
while sum(Cons) ~= m
ce = crit.*x1;
x_i = find(x1);
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
    x1(item(r)) = 0;
end
for j=1:m
         v = A(j,1:end).*x1;
        if (sum(v)<b(j)) && (b(j)>0) 
            Cons(j) = 1;
        elseif (b(j)<0) && (sum(v)>b(j))
            Cons(j) = 1;
        else
            Cons(j) = 0;
        end
end
maxiter = maxiter-1;
if maxiter==0 || sum(x1)<2
    break;
end
end
%Cálculo función objetivo
used_rec = (A*x1');
Z = (W*x1');
x = x1;
%Cálculo recursos usados
if sum(Cons)==m
   fact = 1;
else
    fact = 0;
end

end