
x = [0:2*pi/64:2*pi-2*pi/64];
y = x+1; Lx = 2*pi; Ly = Lx;
[X,Y] = meshgrid(x,y);
u = sin(X)+sin(Y);
[ux, uy] = deri_fft2(u, Lx, Ly);
ux_exact = cos(X);
uy_exact = cos(Y);
figure(1);
subplot(2,2,1); contourf(ux);
subplot(2,2,2); contourf(uy);
subplot(2,2,3); contourf(ux_exact);
subplot(2,2,4); contourf(uy_exact);
[uxx, uxy] = deri_fft2(ux, Lx, Ly);
[uyx, uyy] = deri_fft2(uy, Lx, Ly);
Lap = uxx + uyy;
Lap_exact = -sin(X)-sin(Y);
figure(2);
subplot(2,2,1); contourf(Lap);
subplot(2,2,2); contourf(Lap_exact);
subplot(2,2,2); contourf(Lap_exact);
