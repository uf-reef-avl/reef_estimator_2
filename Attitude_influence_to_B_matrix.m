syms r p y real

B = [1 0 -sin(p);
     0 cos(r) sin(r)*cos(p);
     0 -sin(r) cos(r)*cos(p)]
 
 
A = simplify(inv(A) )
 
