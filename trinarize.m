function data = trinarize(data)

data(data>0) = 1; data(data<0) = -1;