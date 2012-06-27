% Create the diffusion operator
fprintf('Creating diffusion operator...\n');
T=MakeCircleDiffusion(256);   

% Create diffusion wavelets
fprintf('Creating diffusion wavelets...\n');
Tree = DWPTree(T,8,1e-10,struct('Wavelets',false,'ExtendBases',1,'StopFcns',10));

% Create a function to diffuse
F=[zeros(1,floor(size(T,1)/2)),floor(size(T,1)/2):-1:1]'; % Create the vector to take powers of
figure;plot(F);title('Vector to take powers of');

% Take the power 2^p
p = 12;
fprintf('Computing power in compressed fashion...\n');
F_p = DWApplyDyadicPower( Tree, F, p, T );    % Compute a power

% Plot results
figure;subplot(2,1,1);plot(T^(2^(p))*F,'b');hold on;plot(F_p,'r');axis tight;title('Compressed (red) vs. regular (blue) computation');
subplot(2,1,2);plot(log10(abs(T^(2^p)*F-F_p)));axis tight; title('Difference between compressed and regular computation (log scale)');

return;