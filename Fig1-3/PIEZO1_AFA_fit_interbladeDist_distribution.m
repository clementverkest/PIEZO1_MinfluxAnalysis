%% load data
data = result.InterBlades(:,4);

faktor =  size(data,1);

% set limits
MaxSigma1 = 1.5; %1.5 for soma control
MaxSigma2 =3 ; %3 for soma control
MaxSigma3 = 3; %3 for soma control
MaxMean = 40;

%% Options to fit means and/or sigmas
fit_means = true; 
fit_sigmas = true; 
fixed_sigma = 2; % 1.5 for neurite ctl
fixed_means = [15.6, 21.2, 27.1];
% fixed_means = [15.6, 21.2, 27.1];
% Define histogram edges
bin_edges = 5:1:40;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

% Estimate histogram density
[counts, ~] = histcounts(data, bin_edges);
normalized_counts = counts / trapz(bin_centers, counts); % Normalize to make it a probability density


if fit_means && fit_sigmas
    % Initial guesses for means, amplitudes, and standard deviations
    initial_params = [max(normalized_counts), 14, 1, max(normalized_counts), 19, 1, max(normalized_counts), 30, 1];
    
    % Define bounds
    lower_bounds = [0, 0, 0.1, 0, 0, 0.1, 0, 0, 0.1];
    upper_bounds = [Inf, MaxMean, MaxSigma1, Inf, MaxMean, MaxSigma2, Inf, MaxMean, MaxSigma3];
    
    % Define Gaussian sum function with variable means and sigmas
    gaussianSum = @(params, y) gaussian(params(1), params(2), params(3), y) + ...
                                 gaussian(params(4), params(5), params(6), y) + ...
                                 gaussian(params(7), params(8), params(9), y);
elseif fit_means && ~fit_sigmas
    % Initial guesses for means and amplitudes only
    initial_params = [max(normalized_counts), 14, max(normalized_counts), 21, max(normalized_counts), 27];
    
    % Define bounds
    lower_bounds = [0, 0, 0, 0, 0, 0];
    upper_bounds = [Inf, MaxMean, Inf, MaxMean, Inf, MaxMean];
    
    % Define Gaussian sum function with variable means but fixed sigmas
    gaussianSum =@(params, y) gaussian(params(1), params(2), fixed_sigma, y) + ...
                                 gaussian(params(3), params(4), fixed_sigma, y) + ...
                                 gaussian(params(5), params(6), fixed_sigma, y);
elseif ~fit_means && fit_sigmas
    % Fixed means of the Gaussians
    % fixed_means = [14, 21, 27];
    
    % Initial guesses for amplitudes and sigmas only
    initial_params = [max(normalized_counts), 1, max(normalized_counts), 1, max(normalized_counts), 1];
    
    % Define bounds
    lower_bounds = [0, 0.1, 0, 0.1, 0, 0.1];
    upper_bounds = [Inf, MaxSigma1, Inf, MaxSigma2, Inf, MaxSigma3];
    
    % Define Gaussian sum function with fixed means but variable sigmas
    gaussianSum = @(params, y) gaussian(params(1), fixed_means(1), params(2), y) + ...
                                 gaussian(params(3), fixed_means(2), params(4), y) + ...
                                 gaussian(params(5), fixed_means(3), params(6), y);
else
    % Fixed means and sigmas
    % fixed_means = [14, 21, 27];
    
    % Initial guesses for amplitudes only
    initial_params = [max(normalized_counts), max(normalized_counts), max(normalized_counts)];
    
    % Define bounds
    lower_bounds = [0, 0, 0];
    upper_bounds = [Inf, Inf, Inf];
    
    % Define Gaussian sum function with fixed means and sigmas
    gaussianSum = @(params, y) gaussian(params(1), fixed_means(1), fixed_sigma, y) + ...
                                 gaussian(params(2), fixed_means(2), fixed_sigma, y) + ...
                                 gaussian(params(3), fixed_means(3), fixed_sigma, y);
end

% Fit the function using lsqcurvefit with max iterations
opts = optimset('Display', 'off', 'MaxIter', 500); % Suppress output and set max iterations
gaus_params_fit = lsqcurvefit(@(p, x) gaussianSum(p, x), initial_params, bin_centers, normalized_counts, lower_bounds, upper_bounds, opts);

% Extract fitted parameters
if fit_means && fit_sigmas
    fitted_amplitudes = [gaus_params_fit(1), gaus_params_fit(4), gaus_params_fit(7)];
    fitted_means = [gaus_params_fit(2), gaus_params_fit(5), gaus_params_fit(8)];
    fitted_sigmas = [gaus_params_fit(3), gaus_params_fit(6), gaus_params_fit(9)];
elseif fit_means && ~fit_sigmas
    fitted_amplitudes = [gaus_params_fit(1), gaus_params_fit(3), gaus_params_fit(5)];
    fitted_means = [gaus_params_fit(2), gaus_params_fit(4), gaus_params_fit(6)];
    fitted_sigmas = [fixed_sigma, fixed_sigma, fixed_sigma];
elseif ~fit_means && fit_sigmas
    fitted_amplitudes = [gaus_params_fit(1), gaus_params_fit(3), gaus_params_fit(5)];
    fitted_means = fixed_means;
    fitted_sigmas = [gaus_params_fit(2), gaus_params_fit(4), gaus_params_fit(6)];
else
    fitted_amplitudes = [gaus_params_fit(1), gaus_params_fit(2), gaus_params_fit(3)];
    fitted_means = fixed_means;
    fitted_sigmas = [fixed_sigma, fixed_sigma, fixed_sigma];
end

% Calculate proportions of each Gaussian component
total_area = sum(fitted_amplitudes .* fitted_sigmas * sqrt(2 * pi));
proportions = (fitted_amplitudes .* fitted_sigmas * sqrt(2 * pi)) / total_area;

% Display proportions
disp('Proportion of each Gaussian component:');
disp(proportions);

% Display means
disp('Mean of each Gaussian component:');
disp(fitted_means);

% Calculate and display goodness of fit (R-squared)
predicted_counts = gaussianSum(gaus_params_fit, bin_centers);
SS_res = sum((normalized_counts - predicted_counts).^2);
SS_tot = sum((normalized_counts - mean(normalized_counts)).^2);
R_squared = 1 - (SS_res / SS_tot);

disp('Goodness of fit (R-squared):');
disp(R_squared);

% Plot histogram and fitted curve
figure;
bar(bin_centers, faktor*normalized_counts, 'FaceColor', [0.6 0.6 0.6],'FaceAlpha', 0.5, 'EdgeColor', 'k', 'DisplayName', 'Histogram'); hold on;
ax = gca;
ax.TickDir = 'out';
ax.Box = 'off';

% Generate smooth x values for plotting
smooth_x = linspace(min(bin_centers), max(bin_centers), 1000);
smooth_y = gaussianSum(gaus_params_fit, smooth_x);
plot(smooth_x, faktor*smooth_y, 'k', 'LineWidth', 1, 'DisplayName', 'Fitted Gaussians');

% Plot individual Gaussians
for i = 1:3
    plot(smooth_x, gaussian(faktor*fitted_amplitudes(i), fitted_means(i), fitted_sigmas(i), smooth_x), 'Color', [0.8 0 1], 'LineWidth', 1.5, 'DisplayName', ['Gaussian ' num2str(i)]);
end

% legend;
xlabel('interblade distance (nm)'); ylabel('counts');
% title('Fitting Sum of Three Gaussians to Histogram of Data');
% grid on;
