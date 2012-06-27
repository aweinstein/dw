function [vDist,vStats] = TreeDistance( cTree_1, cTree_2, cOpts )

%
% function [vDist,vStats] = TreeDistance( cTree_1, cTree_2, cOpts )
%
% Compute the distance between two trees as a weighted sum of distances between the 
% scaling function subspaces
%
% IN:
%   cTree_1, cTree_2    : the two trees to be compared
%   cOpts               : structure with the following fields:
%           [Factor]    : dilation factor. Default: 2.
%           [Exponent]  : smoothness exponent. Default: 1.
%
% OUT:
%   vDist               : \sum_j Factor^(Exponent*j) dist(cTree_1.ExtBasis,cTree_2.ExtBasis)
%                         where dist is the subspace distance
%   vStats              : structure containing the following fields:
%                           Dist    : vector of length equal to the minimal number of scales in the two trees,
%                                       with the distance between the two scaling function subspaces at each scale
%
% SC:
%   MM :    12/07/08
%
% (c) Copyright Duke University, 2008
%   Mauro Maggioni
%

if nargin<3, cOpts = []; end;
if ~isfield(cOpts,'Factor'),    cOpts.Factor = 2; end;
if ~isfield(cOpts,'Exponent'),  cOpts.Exponent = 1; end;

J = min([size(cTree_1,1),size(cTree_2,1)]);

vStats.Dist(J) = subspace( full(cTree_1{J,1}.ExtBasis), full(cTree_2{J,1}.ExtBasis) );
vDist = cOpts.Factor^(cOpts.Exponent*J)*vStats.Dist(J);

for j = J-1:-1:1,
    if (~isempty(cTree_1{j,2}.ExtBasis)) && (~isempty(cTree_2{j,2}.ExtBasis)),
        % Compute distance between scaling function subspaces at this scale
        vStats.Dist(j) = subspace( full(cTree_1{j,2}.ExtBasis), full(cTree_2{j,2}.ExtBasis) );
    else
        if xor( (~isempty(cTree_1{j,2}.ExtBasis)),(~isempty(cTree_2{j,2}.ExtBasis))),
            vStats.Dist(j) = pi/2;
        else
            vStats.Dist(j) = 0;
        end;
    end;
    % Update the total distance
    vDist = vDist + cOpts.Factor^(cOpts.Exponent*j)*vStats.Dist(j);
end;

return;