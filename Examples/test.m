[Theta Phi Rho] = cart2sph(Points(:,1), Points(:,2), Points(:,3));

%F = cos(Theta.^3);
F = sin(3*Theta.^8) .* exp(-5*Phi.^2);


G = cos(Phi).^3;

idxs = find(-pi/2 <= Phi & Phi <= pi/2 & pi/4<= Theta & Theta <= pi/2);
G(idxs) = 1.4;
