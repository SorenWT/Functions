function [data] = data2ft(matin,srate)
% takes arbitrary data and makes a Fieldtrip data structure out of it
% matin should be either a vector or a matrix in the form n channels x m
% time points
% srate is the sampling rate of the data

if size(matin,2)==1
   matin = horz(matin); 
end

data = struct;
data.trial{1} = matin;
data.time{1} = linspace(0,length(matin)/srate-1/srate,length(matin));
data.fsample = srate;
data.label = cellcat('channel',erase(cellstr(num2str([1:size(data.trial{1},1)]')),' '),'',0)';

