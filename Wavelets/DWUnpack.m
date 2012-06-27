function CoeffList = DWUnpack(Basis, FcnNum)
% function CoeffList = DWUnpack(Basis, [FcnNum])
%
% 
%
%

if ~exist('FcnNum')
   FcnNum = 1;
end

% figure out if we need to fetch more than one function
if length(FcnNum) > 1
end


CoeffList = [];

MaxLevels = size(Basis,1);
MaxNodes  = size(Basis,2);


for Level = 1:MaxLevels
   for Node = 1:MaxNodes
      if ~isempty(Basis{Level, Node})          
         M = size(Basis{Level, Node}, 1);
         if prod(size(Basis{Level,Node}))>1,
            CoeffList = [CoeffList; Level*ones(M,1) Node*ones(M,1) (1:M)' Basis{Level,Node}(:,FcnNum)];
         end;
      end
   end
end

% [junk idxs] = sort(abs(CoeffList(:,4)));
% CoeffList = CoeffList(flipud(idxs), :);
