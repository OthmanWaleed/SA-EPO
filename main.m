% Self-Adaptive Emperor Penguin Optimizer with Multi-Strategy Parameter Adaptation Mechanism for Complex Optimization Problems
% Othman Waleed Khalid, Nor Ashidi Mat Isa , Wei Hong Lim
% School of Electrical and Electronic Engineering, Engineering Campus, Universiti Sains Malaysia
% Faculty of Engineering, Technology and Built Environment, UCSI University, Malaysia
 
% SA-EPO (Benchmarked with CEC2017)
clear; clc; close all;

dimension = 100; % Set the dimension of the problem (10, 30, 50, and 100)

algorithm_names = {
    'SA_EPO_17'
};

results = cell(30, length(algorithm_names)); % Initialize cell array for results
convergence_curves = cell(30, length(algorithm_names)); % Initialize cell array for convergence curves
best_std_stats = cell(30, 2); % Initialize cell array for mean and std dev fitness statistics
Gvalue = (1:30) .* 100;
for test_function_index = 1:30
    opts = struct();
    opts.maxIter = 200; % Maximum number of iterations
    opts.maxFES = 100000; % Maximum number of function evaluations
    opts.dimension = dimension; % Dimension of the problem
    opts.testFunction = test_function_index; % CEC-2017 test function index

    current_best_fitness = zeros(length(algorithm_names), 1);
    current_std_dev_fitness = zeros(length(algorithm_names), 1);

    for i = 1:length(algorithm_names)  
    if test_function_index ~= 2  
        algorithm_name = algorithm_names{i};  
        % Call the algorithm for benchmarking  
        result = feval(algorithm_name, opts);  

        % Store the results  
        results{test_function_index, i} = result;  

        % Store the convergence curve  
        convergence_curves{test_function_index, i} = result.convergenceCurve;  

        % Calculate mean and standard deviation fitness statistics  
        fitness_values = result.convergenceCurve;  
        current_best_fitness(i) = min(fitness_values) - Gvalue(test_function_index);  
        current_std_dev_fitness(i) = std(fitness_values);  

        % Display the results  
        disp(['F ', num2str(test_function_index), ', Best Fit (', algorithm_name, '): ', num2str(current_best_fitness(i))]);  
    end  
end

    % Store mean and standard deviation fitness statistics for the current test function
    best_std_stats{test_function_index, 1} = current_best_fitness;
    best_std_stats{test_function_index, 2} = current_std_dev_fitness;
end

% Prompt to choose the directory for saving the results
result_dir = uigetdir(pwd, 'Select the directory to save the results');

% Export mean and standard deviation fitness data to Excel files in the selected directory
best_fitness_file = fullfile(result_dir, 'best_fitness_data.xlsx');
std_dev_fitness_file = fullfile(result_dir, 'std_dev_fitness_data.xlsx');

xlswrite(best_fitness_file, cell2mat(best_std_stats(:, 1)), 'Sheet1');
xlswrite(std_dev_fitness_file, cell2mat(best_std_stats(:, 2)), 'Sheet1');
% Define unique markers for each algorithm
markers = {'d'};

% Plot the convergence curves
for test_function_index = 1:30
    figure;
    hold on;

    for i = 1:length(algorithm_names)
        algorithm_name = algorithm_names{i};
        mean_val = best_std_stats{test_function_index, 1}(i);
        curve = convergence_curves{test_function_index, i};

        % Plotting every 10th iteration
        x_values = 0:10:(length(curve) - 1) * 10;
        y_values = curve;

        % Determine the color based on the algorithm
        if strcmp(algorithm_name, 'SA_EPO_17')
            color = 'red';
     
        end

        % Plot the curve with the specified color and marker
        plot(x_values, y_values - Gvalue(test_function_index), 'LineWidth', 1, 'Color', color, 'Marker', markers{i});
    end

    hold off;

    xlabel('Iterations', 'FontSize', 16);
    ylabel('Fitness Value', 'FontSize', 16);
    title(['F', num2str(test_function_index)], 'FontSize', 18);
    set(gca, 'FontSize', 16); % Set font size for both x and y axis

    xlim([0, 200]); % Setting x-axis limits
    xticks(0:10:200);

    % Change the legend labels
   % new_legend_labels = {
    %     'SA-EPO'
    %};

    % Place the legend outside the plot
    %legend(new_legend_labels, 'Location', 'eastoutside', 'FontSize', 16);
end
