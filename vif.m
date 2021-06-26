function V = vif(data)

R0 = corrcoef(data,'rows','complete');
V = diag(inv(R0))';