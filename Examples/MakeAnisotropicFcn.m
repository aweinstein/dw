%function F = MakeAnisotropicFcn()

figure; imagesc(T);
figure; plot(lImpedance);

Thetas = linspace(0, 2*pi, 257);
Thetas = Thetas(1:256);


HF = sin(17*Thetas.^5);

CENTER = 196;

Bump1 = zeros(1, 256);
Bump2 = zeros(1, 256);
Spike = zeros(1, 256);

for j=128:256
   Bump1(j)  =  exp ( - 1/850*(j-CENTER)^2);
end

for j=1:128
   Bump2(j)  =  exp ( - 1/850*(j-64)^2);
end


BigBump = zeros(1,256);

for j=1:256
   BigBump(j) = exp( - 1/4000*(j-128)^2);
end

F = sin(2*Thetas.^2).*Bump2+HF.*Bump1 + Spike;

figure; plot(F);
figure; plot(F2);

%[CoeffTree Coeffs] = DiffusionBestBasis(Tree, F');
%[CoeffTree2 Coeffs2] = DiffusionBestBasis(Tree, F2');

%figure;
plot(abs(Coeffs(1:100,1)));
hold on
plot(abs(Coeffs2(1:100,1)),'r');