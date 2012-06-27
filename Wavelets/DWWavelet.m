function Basis = DWWavelet(Coeffs)
% function Basis = DWWavelet(Coeffs)
%
% DWWAVELET extracts the coefficients in the wavelet 


Levels = size(Coeffs,1);
Basis  = cell(size(Coeffs));

level = 1;
while level < Levels & ~isempty(Coeffs{level,1}) 
   Basis{level,2} = Coeffs{level,2};
   level = level+1;
end

Basis(level-1, 1) = Coeffs(level-1,1);
