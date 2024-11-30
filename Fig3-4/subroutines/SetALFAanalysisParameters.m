
function [p]=SetALFAanalysisParameters()

%% set analysis parameters
% set DBSCAN parameters for identifying signals from the same channel
p.MaxDistProtamers = 8;%e-9;
p.MinNumProtamers = 2;

% set DBSCAN parameters for identifying PIEZO channel clusters
p.MaxDistWithinCluster = 65;%;e-9;
p.MinNumChannelsPerCluster = 12;

% set filtering parameters
p.z_Upper_threshold = 7e-07;  % Z-filt threshold
p.z_Lower_threshold = -7e-07;  % Z-filt threshold
p.stdev_trace_threshold = 10; %e-09; % MAX standard deviation
p.loc_per_trace_threshold = 3; % MIN number of localisations per trace
p.efo_threshold = 200000; % MAX efo
p.cfr_threshold = 0.8; % MAX cfr
p.time_threshold = 200000; % MAX recording time

% trimer detection thresholds
p.MaxInterBladeDist=40; %e-9;
p.MinTrimerDist=60; %e-9; % must be greater than p.MaxInterBladeDist
end