function Points = GeneratePoints(Name, NPoints, varargin)
% function Points = GeneratePoints(Name, NPoints, ...)
%
% GENERATEPOINTS generates point clouds on one of a number of special sets,
% including the unit cube and unit sphere.  Along with the name of the graph
% and number of points to generate, GENERATEPOINTS accepts arguments depending
% on the set:
%
%  GeneratePoints('Sphere', NPoints, n)
%     Uniformly distributed points on the sphere S^n) embedded in R^(n+1).
%
%  GeneratePoints('Cube', NPoints, n)
%     Uniformly distributed points on the cube [0,1]^n embedded in R^n.
%
%  GeneratePoints('Disc', NPoints)
%     Uniformly distributed points in unit disc in R^2.
%
%  GeneratePoints('Interval', NPoints)
%     Evenly spaced points on the unit interval [0,1].
%
%  GeneratePoints('Circle', NPoints)
%     Evenly spaced points on the circle S^1, in R^2.
%
%  GeneratePoints('Annulus', NPoints, r1, r2)
%     Uniformly distributed points on the annulus r1 <= r <= r2 in R^2.
%
% In:
%    Name       = string indicating which type of graph to generate
%    NPoints    = number of points
%
% Out:
%    Points     = DxN matrix specifying the N points in R^D which comprise the
%                 generated point cloud
%
% Dependencies:
%    none
%
% Version History:
%    jcb       12/2006        initial version created from MakeExamples.m
%

Points = [];

% shortcut variable
N = NPoints;

% determine the type of graph to generate
if strcmpi(Name, 'Sphere')
   if nargin > 2
      d = varargin{1};
   else
      d = 2;
   end

   Points = randn(d+1, N);
   for j=1:N
      Points(:,j) = Points(:,j)/norm(Points(:,j));
   end

elseif strcmpi(Name, 'Cube')
   if nargin > 2
      d = varargin{1};
   else
      d = 3;
   end

   Points = rand(N, Dim);
elseif strcmpi(Name, 'Disc')

   r = rand(1, N);
   theta = 2*pi*rand(1,N);
   Points = [sqrt(r).*cos(theta); sqrt(r).*sin(theta)];

elseif strcmpi(Name, 'Interval')
   Points = (0:1/N:(1-1/N));
elseif strcmpi(Name, 'Circle')
   Thetas = 0:2*pi/N:2*pi*(1-1/N);
   Points = [cos(Thetas); sin(Thetas)];
elseif strcmpi(Name, 'Annulus')

   if nargin < 4
      fprintf('GeneratePoints.m: missing arguments r1 and r2\b');
      return;
   end

   r1 = varargin{1};
   r2 = varargin{2};

   r = (r2-r1)*rand(1, N)+r1;
   theta = 2*pi*rand(1,N);
   Points = [sqrt(r).*cos(theta); sqrt(r).*sin(theta)];
else
   fprintf('GeneratePoints.m: unknown graph "%s"\n', Name);
   return;
end

