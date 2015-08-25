% https://perswww.kuleuven.be/~u0017946/publications/Papers97/art97a-Saff-Kuijlaars-MI/Saff-Kuijlaars-MathIntel97.pdf

% http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere

N = 7;

s  = 3.6/sqrt(N);
dz = 2.0/(2*N);
long = 0;
z    = 1 - dz/2;

node = NaN(N,3);

for k = 1 : N
    r    = sqrt(1-z*z)
    node(k,:) = [cos(long)*r, sin(long)*r, z];
    z    = z - dz
    long = long + s/r
end

plot3(node(:,1),node(:,2),node(:,3),'.');
