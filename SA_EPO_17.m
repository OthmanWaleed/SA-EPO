% Self-Adaptive Emperor Penguin Optimizer with Multi-Strategy Parameter Adaptation Mechanism for Complex Optimization Problems
% Othman Waleed Khalid, Nor Ashidi Mat Isa , Wei Hong Lim
% School of Electrical and Electronic Engineering, Engineering Campus, Universiti Sains Malaysia
% Faculty of Engineering, Technology and Built Environment, UCSI University, Malaysia
 
% Self-adaptive EPO benchmarked with CEC2017  
function EPO = SA_EPO_17(opts)
    % Initialization
    lb = 0; % Lower bound for optimization
    ub = 1; % Upper bound for optimization

    % Extract parameters from the opts struct
    if isfield(opts, 'maxIter')
        max_Iter = opts.maxIter;
    else
        max_Iter = 200; % Default maximum number of iterations
    end
        
    if isfield(opts, 'maxFES')
       max_FES = opts.maxFES;
    else
       max_FES = 100000; % Default maximum number of iterations
    end

    if isfield(opts, 'dimension')
        dim = opts.dimension;
    else
        dim = 100; % Default dimension
    end

    % Extract the CEC-2017 test function index
    if isfield(opts, 'testFunction')
        func = opts.testFunction;
    else
        error('CEC-2017 test function index not provided.');
    end

    N = 50; % Number of penguins (population size)

    % Initialize variables
    fit = zeros(1, N);
    fitG = inf;
    curve = zeros(1, max_Iter);
    fitcount =0;
    t = 1;

    % The initialization of EPO control parameters
    f_min = 2; % Minimum value of f 
    f_max = 3; % Maximum value of f 
    l_min = 1.5; % Minimum value of l 
    l_max = 2; % Maximum value of l 
    M_min = 1.0; % Minimum value of M
    M_max = 1.50; % Maximum value of M 
    
    if isfield(opts,'T'), max_Iter = opts.T; end
    if isfield(opts,'N'), N = opts.N; end
    if isfield(opts,'M'), M = opts.M; end 
    if isfield(opts,'f'), f = opts.f; end 
    if isfield(opts,'l'), l = opts.l; end 
     
    % Population initialization within the bounds [0, 1]
    X = rand(N, dim); % Random initialization
    X = lb + (ub - lb) * X; % Scale to desired bounds

    % Initialize Strategy Selection Probability (SSP) 
    straNum = 5; % Number of strategies
    SSP = ones(1, straNum) / straNum; % SSP vector for 5 strategies (0.2 each)
    
    LP = 10; % Learning period
    SFlag = zeros(N, straNum); % Success flag matrix
    FFlag = zeros(N, straNum); % Failure flag matrix
    Stotal = [];  % Total success over one learning periods (Vector)
    Ftotal = [];  % Total failure over one learning periods (Vector)
   
    strategy = 1; % Initialize chosen strategy outside the loop
    counts = zeros(1, 5);% Initialize counting the selection of each strategy
 
    SSPN = zeros(50, 5); % Initialize 2D Matrix to store SSP changes for each penguin
    SSPChanges = zeros(max_Iter, N, straNum);% Initialize 3D Matrix to store SSP changes for each iteration
        
    % Initialze countsMatrix to store number of stratgy selection 'counts' at each iteration
    countsMatrix = zeros(max_Iter, straNum);
    
    % Initializing current ciculation and last LP circulations variables
    curcirc=1;
    lastcircLP=0;
    
    % Iterations
    while t <= max_Iter && fitcount <= max_FES
         for i = 1:N
             % Fitness evaluation using the CEC-2017 benchmark function
             fit(i) = cec17_func(X(i, :)', func);
             fitcount=fitcount+1;
         end
         stp=min(fit);
            % Generate radius in [0,1]
            R = rand();
            % Time (7)
            if R > 1
                T0 = 0;
            else
                T0 = 1;
            end
            % Temperature profile (7)
            T = T0 - (max_Iter / (t - max_Iter));
           
         for i = 1:N   
             
             % Update global best and update flags
            if fit(i) < fitG
                fitG = fit(i);  
                Xgb = X(i, :);
                SFlag(i, strategy) = 1; % update success flag
            else
                FFlag(i, strategy) = 1; % update failure flag
            end
            

            % Ckeck LP to update SSP probability 
             if curcirc - lastcircLP == LP
                lastcircLP = curcirc;
                
                % Update total flags for success and failure
                [Stotal, Ftotal, SFlag, FFlag] = totalflags(SFlag, FFlag, straNum, Stotal, Ftotal, N);
                % Calculate SSP for 5 strategies based on flags of success and failure
                [SSP, Stotal, Ftotal] = updateSSP(Stotal, Ftotal, SSP, straNum);
                % Select current strategy by using roulette wheel selection
                [strategy] = rouletteWheelSelection(SSP, straNum); 
             end
           
            curcirc = curcirc + 1; % increament circulation for LP tracking
            SSPN (i, :) = SSP; % Update probabilities matrix inside N loop 
            fitness = fit(i);
            
            % function to adaptively change control parameters using selected strategy
            % Also counts how many times each strategy has been selected
            [f, l, M, counts] = updateParameters(strategy, f_min, f_max, l_min, l_max, M_min, M_max, fitness, counts);
             
            for d = 1:dim
                % Pgrid 
                P_grid = abs(Xgb(d) - X(i, d));
                % Vector A 
                A = (M * (T + P_grid) * rand()) - T;
                % Vector C 
                C = rand(); 
                % Compute function S 
                S = sqrt(f * exp(t / l) - exp(-t))^2;
                % Distance 
                Dep = abs(S * Xgb(d) - C * X(i, d));
                % Position update 
                X(i, d) = Xgb(d) - A * Dep;  
             end
             
            % Boundary
            XB = X(i, :);
            XB(XB > ub) = ub;
            XB(XB < lb) = lb;
            X(i, :) = XB;
           
        end
           
        % The convergence curve
          if t==1
            curve(t) =stp;
        else
        curve(t) = fitG;
          end
          
       
        SSPChanges(t, :, :) = SSPN; % 3D matrix to capture SSP changes at every iteration
        SSPN = zeros(N, straNum); % reset storage of SSP changes in N loop
         
         % Store 'counts' in the matrix for each iteration
         countsMatrix(t, :) = counts;        
         counts = zeros(1, 5);% reset counting storage of all strategies
           
          if fitcount>=max_FES
             break;
          end
         t = t + 1;
    end
    
% -----------------------------------------------------------------------------%
% Plot the probability of each strategy over iterations 'SSPChanges'
%------------------------------------------------------------------------------%

% Store results
    EPO.bestFitness = fitG;
    EPO.convergenceCurve = curve;
end

%%%-------------- helper functions -----------------%%%%

function [strategy] = rouletteWheelSelection(SSP, num)
 
    straNum = length(SSP);  
    Select = zeros(1, num); % Array to store selected strategies
    r = rand(1, num); % Array of random values between 0 and 1

    % Loop to perform roulette wheel selection
    for i = 1:num
        sumSSP = 0;  
        r2 = randperm(straNum); % Generate a random permutation of strategy indices
        j = r2(1); % Initialize the index for strategy selection

        while sumSSP < r(i)
            sumSSP = sumSSP + SSP(mod(j - 1, straNum) + 1); % Accumulate strategy selection probabilities
            j = j + 1; % Move to the next strategy index
        end

        % Record the selected strategy index
        Select(i) = mod(j - 2, straNum) + 1;
    end

    tempSelect = zeros(1, straNum);

     % Count the occurrences of each strategy in the selection
    for i = 1:straNum
        tempSelect(i) = sum(Select == i);
    end

    % Select the strategy with the highest count
    [~, sortedIndices] = sort(tempSelect, 'descend');
    strategy = sortedIndices(1);
 
end

 
function [f, l, M, counts] = updateParameters(strategy, f_min, f_max, l_min, l_max, M_min, M_max, fitness, counts)
    % Validate strategy
    if ~isscalar(strategy) || ~isnumeric(strategy) || strategy < 1 || strategy > 5
        error('Strategy must be an integer between 1 and 5.');
    end
 
    % Calculate min and max fitness values
    minFit = min(fitness);
    maxFit = max(fitness);

    % Update parameters based on the selected strategy
    switch strategy
        case 1 % Linear
            f = f_min + (f_max - f_min) * (1 - (fitness - minFit) / (maxFit - minFit));
            l = l_min + (l_max - l_min) * (1 - (fitness - minFit) / (maxFit - minFit));
            M = M_min + (M_max - M_min) * (1 - (fitness - minFit) / (maxFit - minFit));
        case 2 % Quadratic
            f = f_min + (f_max - f_min) * (1 - (fitness - minFit) / (maxFit - minFit)).^2;
            l = l_min + (l_max - l_min) * (1 - (fitness - minFit) / (maxFit - minFit)).^2;
            M = M_min + (M_max - M_min) * (1 - (fitness - minFit) / (maxFit - minFit)).^2;
        case 3 % Exponential
            f = f_min + (f_max - f_min) * exp(-(fitness - minFit) / (maxFit - minFit));
            l = l_min + (l_max - l_min) * exp(-(fitness - minFit) / (maxFit - minFit));
            M = M_min + (M_max - M_min) * exp(-(fitness - minFit) / (maxFit - minFit));
        case 4 % Logarithmic
            f = f_min + (f_max - f_min) * log((fitness - minFit) / (maxFit - minFit) + 1);
            l = l_min + (l_max - l_min) * log((fitness - minFit) / (maxFit - minFit) + 1);
            M = M_min + (M_max - M_min) * log((fitness - minFit) / (maxFit - minFit) + 1);
        case 5 % Trigonometric (Sine)
            f = f_min + (f_max - f_min) * sin(pi/2 * (fitness - minFit) / (maxFit - minFit));
            l = l_min + (l_max - l_min) * sin(pi/2 * (fitness - minFit) / (maxFit - minFit));
            M = M_min + (M_max - M_min) * sin(pi/2 * (fitness - minFit) / (maxFit - minFit));
    end
    
     % Update counts of each strategy for plotting later
     counts(strategy) = counts(strategy) + 1;
    
    % Ensure parameters stay within defined bounds
    f = max(min(f, f_max), f_min);
    l = max(min(l, l_max), l_min);
    M = max(min(M, M_max), M_min);
end

function [Stotal,Ftotal,SFlag,FFlag]=totalflags(SFlag,FFlag,straNum,Stotal,Ftotal,N)
TempStotal=sum(SFlag,1);
Stotal=[Stotal;TempStotal];
TempFtotal=sum(FFlag,1);
Ftotal=[Ftotal;TempFtotal];
SFlag=zeros(N,straNum);
FFlag=zeros(N,straNum);
end

function [SSP, Stotal, Ftotal] = updateSSP(Stotal, Ftotal, ~,straNum)

     tempSum1=sum(Stotal,1);    
        tempSum2=tempSum1;  
        [~,j]=size(tempSum1);  
        for k=1:j
           if (tempSum1(k)==0)  
              tempSum2(k)=1;  
           end
        end
        if (sum(tempSum1(1,:),2)==0)  
        SSP = rand(1, straNum);
        SSP = SSP / sum(SSP);  
        else
        SSP=tempSum1./(tempSum2+sum(Ftotal,1));
        SSP=SSP./sum(SSP,2);
        end
        Stotal=[]; % Reset the total success matrix
        Ftotal=[]; % Reset the total failure matrix
end