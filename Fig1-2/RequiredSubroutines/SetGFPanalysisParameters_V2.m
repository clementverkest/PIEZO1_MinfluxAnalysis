
function [p]=SetGFPanalysisParameters_V2()

%% set analysis parameters
% set DBSCAN parameters for identifying signals from the same channel
p.MaxDistProtamers = 25;%e-9;
p.MinNumProtamers = 2;

% set DBSCAN parameters for identifying PIEZO channel clusters
p.MaxDistWithinCluster = 135;%;e-9;
p.MinNumChannelsPerCluster = 12;

% set filtering parameters
p.z_Upper_threshold = 7e-07;  % Z-filt threshold
p.z_Lower_threshold = -7e-07;  % Z-filt threshold
p.stdev_trace_threshold = 10; %e-09; % MAX standard deviation
p.loc_per_trace_threshold = 3; % MIN number of localisations per trace
p.efo_threshold = 200000; % MAX efo
p.cfr_threshold = 0.5; % MAX cfr
p.time_threshold = 2000000; % MAX recording time

end