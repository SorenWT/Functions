function tabout = ttimes(tabin,mult)

tabout = array2table(tabin{:,:}.*mult,'VariableNames',tabin.Properties.VariableNames);