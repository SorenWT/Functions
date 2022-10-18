function [ T,P,pcvar ] = nipals(X,a,it,tol )
%Nipals algorithm for Principle Component Analysis
%
% Input:
%   X: a matrix with rows obseravations and columns variables. 
%   a: number of components desired in the output. 
%   it: maximum number of iterations allowed to find one component. default 1000.
%   tol: tolerance value to assese convergence. default 1e-4
%
% Output:
%   T: scores or transformed coordinates.
%   P: loadings or coefficients.
%   pcvar: variance captured in each output component.
%
% If X is m*n matrix, then T is m*a, P is n*a
%
% This function is written largely based on nipals function from R chemometrics
% package. It was frustrating that Matlab did not has a robust nipals
% function so I wrote one for my project and share here.
%
% Author: Qiaonan Duan, 6/7/2013, MSSM.
%
% Updates:
% Line 48 (6/10/2013): prec = thnew-th; --> prec = (thnew-th)'*(thnew-th);
%
if nargin == 2
		it = 1000;
		tol = 1e-4;
	elseif nargin == 3
		tol = 1e-4;
	end
	[obsCount,varCount] = size(X);
	Xh = X - repmat( mean(X,1), obsCount, 1 );
	T = zeros(obsCount,a);
	P = zeros(varCount,a);
	pcvar = zeros(varCount,1);
	varTotal = sum(var(Xh));
	currVar = varTotal;
	nr = 0;
	for h = 1:a
		th = Xh(:,1);
		ende = false;
		while(~ende)
			nr = nr+1;
			ph = Xh'*th/(th'*th);
			ph = ph/norm(ph);
			thnew = Xh*ph/(ph'*ph);
			prec = (thnew-th)'*(thnew-th);
			th = thnew;
			if prec <= tol^2
				ende = true;
			elseif it <= nr
				ende = true;
				disp('Iteration stops without convergence')
			end
		end
		Xh = Xh-th*ph';
		T(:,h) = th;
		P(:,h) = ph;
		oldVar = currVar;
		currVar = sum(var(Xh));
		pcvar(h) = ( oldVar - currVar )/varTotal;
		nr = 0;
end
