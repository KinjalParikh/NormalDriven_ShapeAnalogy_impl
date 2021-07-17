[Vi,Fi] = readOBJ('spot.obj');
Ni = vertexNormal(triangulation(Fi,Vi));    %returns unit normal vector

[Vs,Fs] = subdivided_sphere(3);
Ns = vertexNormal(triangulation(Fs,Vs));

[Vt, Ft] = readOBJ('cube.obj');
Nt = vertexNormal(triangulation(Ft,Vt));

[Mt , It] = min(acos((Ns*Nt')), [], 2);
Nst = Nt(It, :);    % Fs -> Ft

[Mi , Ii] = min(acos((Ni*Ns')), [], 2);
Nit = Nst(Ii, :);   % Fi -> Fs -> Ft     Paper: T


V = Vi;

%ARAP Pre-computation
Q = cotmatrix(V, Fi);
preF = [];
[~,K] = arap_rhs(V,Fi,[],'Energy','spokes-and-rims');


%Computing matrix E, W, A
V_F = cell(size(Vi, 1));         %all the triangles in the ring in the form of vertex list
for f = 1:size(Fi, 1)
    for i = 1:3
        V_F{Fi(f, i)} = [V_F{Fi(f, i)}; Fi(f, 1), Fi(f, 2), Fi(f, 3), f];
    end
end
A = massmatrix(Vi, Fi, 'voronoi');

lambda = 500;
M1M2 = cell(size(Vi, 1));
C = cotangent(Vi,Fi);
for i = 1:size(Vi, 1)
   v1 = V_F{i}(:,1); v2 = V_F{i}(:,2); v3 = V_F{i}(:,3); 
   s12 = Vi(v2,:) - Vi(v1,:);
   s13 = Vi(v3,:) - Vi(v1,:);
   s23 = Vi(v3,:) - Vi(v2,:);
   E = [s23' s13' s12'];
   
   w  = C(V_F{i}(:,4), :);
   W = diag(reshape(w, size(w, 1)*3, 1));
   M1 = [E Ni(i, :)'];
   M2 = [W, zeros(size(W, 1), 1); [zeros(1, size(W, 2)), lambda*A(i, i)]];
   M1M2{i} = M1*M2;    
end

E_ = cell(size(Vi, 1));
Rall = zeros(3, 3, size(Vi, 1));
b = 1000;
bc = Vi(b, :);
for i = 1:100
    %local
    for k = 1:size(Vi, 1)
       v1 = V_F{k}(:,1); v2 = V_F{k}(:,2); v3 = V_F{k}(:,3); 
       s12 = V(v2,:) - V(v1,:);
       s13 = V(v3,:) - V(v1,:);
       s23 = V(v3,:) - V(v2,:);
       E_{k} = [s23; s13; s12];
    end

    for k = 1:size(Vi, 1)
       M3 = [E_{k}; Nit(k, :)];
       X = M1M2{k}*M3;
       [U_svd,S,V_svd] = svd(X);
       r_k = V_svd*U_svd';
       if det(r_k)<0
           r_k = -r_k;
       end    
       Rall(:, :, k) = r_k;
    end
    
    %global
    R = reshape(permute(Rall,[3 1 2]),size(Vi, 1)*3*3, 1); 
    rhs = reshape(K * R,[size(K, 1)/3 3]);
    VPre = V;
    [V, preF] = min_quad_with_fixed(Q,rhs,b,bc,[],[],preF);

    % stopping criteria
    dVpre = sqrt(sum((V - VPre).^2,2));
    dVi = sqrt(sum((V - Vi).^2,2));
    reldV = max(dVpre) / max(dVi);
    if reldV <  1e-3
        break;
    end
end

t = tsurf(Fi,V, 'EdgeColor', 'black');shading interp; axis equal;