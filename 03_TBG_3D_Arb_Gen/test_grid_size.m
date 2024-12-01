Lx = 1e-3;     dx = 0.05e-6;
Ly = 1e-3;     dy = 0.05e-6;
x = (-.5*Lx: dx :.5*Lx);
y = (-.5*Ly: dy :.5*Ly);

[xx, yy] = meshgrid(x, y);