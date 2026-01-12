function [bigtblout] = prepare_tbl_mlm(tblin,clustervar)

tbl_cgm = nancenter(tblin);
[tbl_cwc,tbl_grpmean] = nancenter(tblin,clustervar);
tbl_grpmean = nancenter(tbl_grpmean);

tbl_cgm.Properties.VariableNames = strcat(tbl_cgm.Properties.VariableNames,'_cgm');
tbl_cwc.Properties.VariableNames = strcat(tbl_cwc.Properties.VariableNames,'_cwc');

bigtblout = [tblin tbl_cgm tbl_cwc tbl_grpmean];

