function [weightsout,comps] = rot_sparse(data,weights,oblique)

if nargin < 3
    oblique = false; 
end

tmp = triu(ones(size(weights,2))-eye(size(weights,2)));
ntheta = sum(sum(tmp));

thetas = ga(@(theta)minfun(data,weights,theta),ntheta,[],[],[],[],zeros(ntheta,1),2*pi*ones(ntheta,1),...
    [],optimoptions('ga','Display','iter','MaxGenerations',30));

[a,b,weightsout] = minfun(data,weights,thetas);

comps = data*weightsout;

if oblique
    target = sign(weightsout) .* abs(weightsout).^4;
    
    [weightsout, T] = procrustes(weights, target, 'oblique');
    comps = data*weightsout;
end

end


function [varoptim,T,B] = minfun(data,A,thetas)

tmp = triu(ones(size(A,2))-eye(size(A,2)));

thetas_mat = tmp;
thetas_mat(find(thetas_mat)) = thetas;

[B,T] = dorot(A,eye(size(A,2)),thetas_mat);

comps = data*B; comps = comps./vecnorm(comps')';
varoptim = -nansum(var(comps.^2,[],2));

end

function [C,Ceq] = constrfun(data,A,thetas)

tmp = triu(ones(size(A,2))-eye(size(A,2)));

thetas_mat = tmp;
thetas_mat(find(thetas_mat)) = thetas;

[B,T] = dorot(A,eye(size(A,2)),thetas_mat);

B = B.*sign(B(1,:));

if sum(sum(B>=0))>0.9
   Ceq = zeros(size(thetas));
else
    Ceq = ones(size(thetas));
end
C = -ones(size(thetas));

end


function [B,T] = dorot(B,T,theta)

m = size(B,2);

for i = 1:(m-1)
    for j = (i+1):m
        Bi = B(:,i);
        Bj = B(:,j);
        Tij = [cos(theta(i,j)) -sin(theta(i,j)); sin(theta(i,j)) cos(theta(i,j))];
        B(:,[i,j]) = B(:,[i,j]) * Tij;
        T(:,[i,j]) = T(:,[i,j]) * Tij;
    end
end

end


function [B, T] = procrustes(A, target, type)
%PROCRUSTES Procrustes rotation of FA or PCA loadings.
[d, m] = size(A);

if nargin < 2 || isempty(target)
    error(message('stats:rotatefactors:TargetRequired', 'procrustes'));
elseif any(size(target) ~= [d m])
    error(message('stats:rotatefactors:InputSizeMismatch'));
end
if nargin < 3 || isempty(type)
    type = 'oblique';
else
    typeNames = {'oblique','orthogonal'};
    type = internal.stats.getParamVal(type,typeNames,'Type');
end

% Orthogonal rotation to target
switch type
case 'orthogonal'
    [L, ~, M] = svd(target' * A);
    T = M * L';

% Oblique rotation to target
case 'oblique'
    % LS, then normalize
    T = A \ target;
    T = T * diag(sqrt(diag((T'*T)\eye(m)))); % normalize inv(T)
end
B = A * T;

end
