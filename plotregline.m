function f = plotregline(a,b,varargin)

argsin = varargin;

B = regress(b,[ones(size(a)) a]);

f = plot(linspace(min(a),max(a),1000),B(1)+B(2)*linspace(min(a),max(a),1000),argsin{:});
