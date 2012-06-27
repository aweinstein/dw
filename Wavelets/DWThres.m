function Basis = DWThres( Basis, Opts )

for j = 1:size(Basis,1),
    for k = 1:size(Basis,2),
        Basis{j,k} = Basis{j,k}.*(abs(Basis{j,k})>Opts.Thres);
    end;
end;

return;