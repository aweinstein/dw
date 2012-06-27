function F = DWRecon(Tree, Basis)
% function F = DWRecon(Tree, Basis)
%
%

F = [];


Levels    = size(Tree,1);
MaxNodes  = size(Tree,2);

for j=Levels-1:-1:1
    for r=1:MaxNodes

      LeftNode  = 2*r-1;
      RightNode = 2*r;

      if LeftNode <= MaxNodes & RightNode <= MaxNodes
         if ~isempty(Basis{j+1,LeftNode})
            Basis{j, r} = Tree{j+1,LeftNode}.Basis*Basis{j+1,LeftNode};

            if ~isempty(Basis{j+1, RightNode})
               Basis{j,r} = Basis{j,r}+Tree{j+1, RightNode}.Basis*Basis{j+1, RightNode};
            end
         else
             if ~isempty(Basis{j+1, RightNode})
               Basis{j,r} = Tree{j+1, RightNode}.Basis*Basis{j+1, RightNode};
             end
         end

      end
    end
end


if ~isempty(Basis{1,1})

    F = Tree{1,1}.Basis*Basis{1,1};
    if ~isempty(Basis{1,2})
        F = F + Tree{1,2}.Basis*Basis{1,2};
    end
else
    F = Tree{1,2}.Basis*Basis{1,2};
end


%F = Tree{1,1}.Basis*Basis{1,1}' ^+ Tree{2,1}.Basis'*Basis{1,2};

%F = Coeffs{1,1};
%F = Tree{1,1}.Basis*Coeffs{1,1};