function Q = gs(A)

TIME = cputime;

Q = zeros(size(A));
N = size(A,2);

Projections = 5;

for j=1:N
   fprintf('%5d/%5d', j,N);
   Q(:,j) = A(:,j);
   for k=1:Projections
      ips = Q(:,1:j-1)'*Q(:,j);
      Q(:,j) = Q(:,j)-Q(:,1:j-1)*ips;
   end
   Q(:,j) = Q(:,j) / norm(Q(:,j));
   fprintf('\b\b\b\b\b\b\b\b\b\b\b');
end

fprintf('%g seconds\n', cputime-TIME);