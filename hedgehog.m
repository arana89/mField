function [xx, yy, zz, Hx, Hy, Hz] = hedgehog(n)
% simulate 3-D vector hedgehog of array size n

vec = linspace(-1,1,n);
vec = 1:n;
nc = round((n+1)/2);
vec = vec - nc;

[xx,yy,zz] = meshgrid(vec,vec,vec);

% [th,phi,~] = cart2sph(xx,yy,zz);
% 
% Hx = sin(th).*cos(phi);
% Hy = sin(th).*sin(phi);
% Hz = cos(th);
% 
Hx = xx;
Hy = yy;
Hz = zz;

H_norm = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

Hx = Hx ./ H_norm;
Hy = Hy ./ H_norm;
Hz = Hz ./ H_norm;

Hx(isnan(Hx)) = 0;
Hy(isnan(Hy)) = 0;
Hz(isnan(Hz)) = 0;

% figure;quiver3(xx,yy,zz,Hx,Hy,Hz); axis equal;

end