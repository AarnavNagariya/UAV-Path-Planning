function runMultiUAVPlanner(numUAVs)
% RUNMULTIUAVPLANNER Runs path planning for multiple UAVs in 3D space
%   numUAVs - Number of UAVs to plan paths for (default: 3)
%
% This function orchestrates path planning for multiple UAVs by setting different
% start and end points for each UAV and calling the runUAVABC4 function.
%
% Example usage:
%   runMultiUAVPlanner(3)     % Plan paths for 3 UAVs
%   runMultiUAVPlanner()      % Use default of 3 UAVs

global boundary setstart setfinal node delta_H;

% Set default number of UAVs if not specified
if nargin < 1
    numUAVs = 3;
end

% Create a figure with a title for the simulation
figure('Name', 'Multi-UAV 3D Path Planning', 'NumberTitle', 'off', ...
       'Position', [50, 50, 900, 700]);
hold on;

% Set common parameters for all UAVs
boundary = [500, 0];    % Boundary of the environment [max, min]
node = 4;               % Number of intermediate waypoints
delta_H = [50, 50, 50]; % Altitude offsets for each UAV

% Define start and end points for each UAV
% These are sample positions - adjust as needed for your specific scenario
uavStartPoints = [
    100, 100, 100;    % UAV 1 start
    150, 50, 120;     % UAV 2 start
    80, 140, 80;      % UAV 3 start
    200, 80, 110;     % UAV 4 start
    120, 200, 90;     % UAV 5 start
];

uavEndPoints = [
    400, 400, 150;    % UAV 1 end
    350, 450, 140;    % UAV 2 end
    420, 380, 160;    % UAV 3 end
    380, 320, 170;    % UAV 4 end
    450, 350, 180;    % UAV 5 end
];

% Ensure we have enough start/end points for the requested number of UAVs
if numUAVs > size(uavStartPoints, 1)
    error('Not enough predefined start/end points for %d UAVs. Maximum is %d.', ...
          numUAVs, size(uavStartPoints, 1));
end

% Create obstacle positions (shared among all UAVs)
% Format: center coordinates (x,y,z) and radius
danger_xi = [180, 250, 320];
danger_yi = [200, 250, 150];
danger_zi = [140, 180, 120];
danger_ri = [50, 60, 45];

% Make these available globally
global danger_xi danger_yi danger_zi danger_ri;

% Create a colormap for UAV paths
uavColors = lines(numUAVs);

% Initialize arrays to store all paths for later collision checking
allPaths = cell(1, numUAVs);

% Plan paths for each UAV sequentially
for i = 1:numUAVs
    % Set start and end points for this UAV
    setstart = uavStartPoints(i, :);
    setfinal = uavEndPoints(i, :);
    
    % Display progress
    fprintf('\n===== Planning path for UAV #%d =====\n', i);
    fprintf('Start: [%.1f, %.1f, %.1f]\n', setstart(1), setstart(2), setstart(3));
    fprintf('End: [%.1f, %.1f, %.1f]\n', setfinal(1), setfinal(2), setfinal(3));
    
    % Call the path planner for this UAV
    path = runUAVABC4(i);
    
    % Store the path for collision checking
    allPaths{i} = path;
end

% Check for potential collisions between UAV paths
checkPathCollisions(allPaths, 20); % 20 is the minimum safety distance between UAVs

% Create a legend for the UAVs
legendLabels = cell(1, numUAVs);
for i = 1:numUAVs
    legendLabels{i} = sprintf('UAV #%d', i);
end
legend(legendLabels, 'Location', 'northeast');

% Add title and finalize the figure
title('3D Multi-UAV Path Planning');
grid on;
view(45, 30);
axis equal;

% Keep the figure open
hold on;

end

function checkPathCollisions(paths, safetyDistance)
% Check for potential collisions between UAV paths
% paths: Cell array of UAV paths
% safetyDistance: Minimum safe distance between UAVs

numUAVs = length(paths);
collisionFound = false;

for i = 1:numUAVs
    for j = (i+1):numUAVs
        path1 = paths{i};
        path2 = paths{j};
        
        % Sample points along both paths with higher resolution for better detection
        sampledPath1 = samplePath(path1, 50);
        sampledPath2 = samplePath(path2, 50);
        
        % Check distances between all points in the sampled paths
        minDistance = inf;
        for p1 = 1:size(sampledPath1, 1)
            for p2 = 1:size(sampledPath2, 1)
                dist = norm(sampledPath1(p1, :) - sampledPath2(p2, :));
                if dist < minDistance
                    minDistance = dist;
                    point1 = sampledPath1(p1, :);
                    point2 = sampledPath2(p2, :);
                end
            end
        end
        
        % Report if paths come too close to each other
        if minDistance < safetyDistance
            collisionFound = true;
            fprintf('\nWARNING: Potential collision between UAV #%d and UAV #%d\n', i, j);
            fprintf('Minimum distance: %.2f (safety threshold: %.2f)\n', minDistance, safetyDistance);
            
            % Visualize collision point with a red sphere
            [x, y, z] = sphere(20);
            collision_point = (point1 + point2) / 2;
            surf(5*x + collision_point(1), 5*y + collision_point(2), 5*z + collision_point(3), ...
                'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            text(collision_point(1), collision_point(2), collision_point(3) + 10, ...
                'COLLISION RISK', 'Color', 'r', 'FontWeight', 'bold');
        end
    end
end

if ~collisionFound
    fprintf('\nNo potential collisions detected. All paths maintain at least %.2f distance.\n', safetyDistance);
end

end

function sampledPath = samplePath(path, numSamples)
% Create a more densely sampled version of the path for better collision detection

% If path has only a few points, use linear interpolation
if size(path, 1) < 3
    sampledPath = path;
    return;
end

% Generate parameter values for the path
t = linspace(0, 1, size(path, 1));
t_new = linspace(0, 1, numSamples);

% Interpolate each coordinate
sampledPath = zeros(numSamples, 3);
for dim = 1:3
    sampledPath(:, dim) = interp1(t, path(:, dim), t_new, 'spline');
end

end