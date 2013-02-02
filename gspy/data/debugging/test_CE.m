xmin = 0; ymin = 0; xmax = 4; ymax =4; nx = 25; ny = 25;
fspace = fundefn('spli', [ nx ny], [xmin ymin], [xmax ymax]);

x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
z = zeros(1, nx * ny);

n = 0;
for xx=1:nx
    for yy=1:ny
        n = n + 1;
        domain(n, :) = [x(xx) y(yy)];
        z(n) = sin(x(xx)) - cos(y(yy) .^ 2);
    end
end

c = funfitxy(fspace, domain, z');

X = linspace(xmin, xmax, 100);
Y = linspace(ymin, ymax, 100);

ni =0;

for ix=1:100
    for iy=1:100
        ni = ni + 1;
        Z(ni) = funeval(c,fspace,[X(ix) Y(iy)]);
        Zreal(ni) = sin(X(ix)) - cos(Y(iy) .^ 2);
    end
end


max_abs_err = max(abs(Zreal - Z));
disp(max_abs_err);
