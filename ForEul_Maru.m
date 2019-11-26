function [sol,y,te,ye,ie] = ForEul_Maru(func,tspan,y0,options)
% set the step size for this using options.InitialStep (since the size
% never changes)

if ~isempty(options.InitialStep)
    stepsize = options.InitialStep;
else
    stepsize = 1;
end

if ~isempty(options.Sigma)
    sigma = options.Sigma;
else
    sigma = 1;
end

t = linspace(tspan(1),tspan(2),diff(tspan)/stepsize)';
y(:,1) = y0;
te = [];
ye = [];
ie = [];

deriv = 0;

for i = 1:(length(t)-1)
    if ~isempty(options.Events)
        [pos,isterm,dirn] = options.Events(t(i),y(i,:));
        if (pos == 0)
            te = [te t(i)];
            ye = [ye y(i,:)];
            ie = [ie i];
            if (deriv*dirn>=0) && isterm
                t = t(1:i);
                if nargout == 1
                    sol.t = t;
                    sol.y = y;
                    sol.te = te;
                    sol.ye = ye;
                    sol.ie = ie;
                else
                    sol = t;
                end
                return
            end
        end
    end
    
    
    deriv = func(t(i),y(:,i));
    y(:,i+1) = y(:,i) + stepsize*(deriv+(sigma.*randn(size(deriv))./sqrt(stepsize)));
    
end


if nargout == 1
    sol.t = t;
    sol.y = y;
    sol.te = te;
    sol.ye = ye;
    sol.ie = ie;
else
    sol = t;
end
