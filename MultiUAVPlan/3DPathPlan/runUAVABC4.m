function POS=runUAVABC4(obj)
global boundary setstart setfinal node delta_H setstart_2 setfinal_2 danger_xi danger_yi danger_zi danger_ri;
maxCycle=20;% Maximum algorithm iterations
%/* Problem specific variables*/

% Check if obj is provided, if not, default to 1
if nargin < 1 || isempty(obj)
    obj = 1;
end

% Check if figure exists and create or use it
if ~ishold
    figure;
    hold on;

    % Visualize 3D terrain environment
    try
        SETenvironment;
        surf(X,Y,Z);
        box on;
        rotate3d on;
        xi=linspace(0,500,100);
        yi=linspace(0,500,100);
        [XI,YI]=meshgrid(xi,yi);
        ZI=interp2(X,Y,Z,XI,YI,'cubic');
        surf(XI,YI,ZI); % Smoothed surface + contours
        colormap('summer'); % Using a nice colormap for terrain
        alpha(0.7); % Make terrain semi-transparent
        
        % Set proper view and lighting
        view(45, 30);
        lighting phong;
        grid on;
    catch
        warning('Could not visualize terrain environment. SETenvironment function might be missing.');
    end

    % Visualize obstacles if they exist
    if exist('danger_xi', 'var') && ~isempty(danger_xi) && ...
       exist('danger_yi', 'var') && ~isempty(danger_yi) && ...
       exist('danger_zi', 'var') && ~isempty(danger_zi) && ...
       exist('danger_ri', 'var') && ~isempty(danger_ri)
        try
            [x,y,z]=sphere(40);
            for k=1:size(danger_xi,2)
                surf(danger_ri(k)*x+danger_xi(k), danger_ri(k)*y+danger_yi(k), danger_ri(k)*z+danger_zi(k), 'FaceColor', [0.7 0.1 0.1]);
                alpha(0.6); % Make obstacles semi-transparent
            end
        catch
            warning('Could not visualize obstacles.');
        end
    end
    
    % Set up the figure properties
    title('3D Multi-UAV Path Planning');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Z Coordinate (Altitude)');
    axis([boundary(2) boundary(1) boundary(2) boundary(1) 0 500]);
end

D=3; 
NP=40; 
FoodNumber=NP/2; 
limit=NP*D;
ub=ones(1,D)*(boundary(1)-10); 
lb=ones(1,D)*(boundary(2)+10);
setfinal_2=setfinal+[0 0 delta_H(1)];
setstart_2=setstart+[0 0 delta_H(1)];
runtime=node+1;
objfun='UAV2'; 

GlobalMins=zeros(1,runtime);
GlobalParams_s=[setstart_2;zeros(node+1,3);setfinal_2];

for r=1:2*runtime
Range = repmat((ub-lb),[FoodNumber 1]);
Lower = repmat(lb, [FoodNumber 1]);
Foods = rand(FoodNumber,D) .* Range + Lower;

%====================================================================================%
if r<=node/2 % Step1: Plan odd position nodes in pairs, distributing by X coordinates
   X_ave=setstart_2(1):((setfinal_2(1)-setstart_2(1))/(node/2+1)):setfinal_2(1);
   X_ave=X_ave(2:end-1);
   if rem(r,2)~=0
     Foods(:,1)=X_ave((r+1)/2);
     ObjVal=feval(objfun,Foods,GlobalParams_s(r,:),GlobalParams_s(node+4-r,:));   
   else
     Foods(:,1)=X_ave(node/2+1-r/2);    
     ObjVal=feval(objfun,Foods,GlobalParams_s(r+1,:),GlobalParams_s(node+5-r,:));   
   end
else
    if r>node/2&&r<=runtime % Step2: Plan even position nodes individually, distributing by Y coordinates based on nodes from Step1
   feckParams_s=[0;0;GlobalParams_s(1:(end-2),2)];
   feckParams_s=(GlobalParams_s(:,2)-feckParams_s)/2;
   feckParams_s=[feckParams_s(3:end);0;0];
   Y_ave=GlobalParams_s(:,2)+feckParams_s;
   Foods(:,2)=Y_ave(2*(r-node/2)-1);
   ObjVal=feval(objfun,Foods,GlobalParams_s(2*(r-node/2)-1,:),GlobalParams_s(2*(r-node/2)+1,:));   
    elseif r>runtime&&r<=runtime+node/2 % Step3: Refine odd position nodes individually, distributing by X coordinates based on nodes from Step2
   feckParams_s=[0;GlobalParams_s(4:end,1);0;0];     
   feckParams_s=(feckParams_s-GlobalParams_s(:,1))/2;  
   X_ave=GlobalParams_s(:,1)+feckParams_s;
   Foods(:,1)=X_ave(2*(r-runtime));
   ObjVal=feval(objfun,Foods,GlobalParams_s(2*(r-runtime),:),GlobalParams_s(2*(r-runtime)+2,:));   
    else % Step4: Refine even position nodes individually, distributing by Y coordinates based on nodes from Step3
   feckParams_s=[0;0;GlobalParams_s(1:(end-2),2)];
   feckParams_s=(GlobalParams_s(:,2)-feckParams_s)/2;
   feckParams_s=[feckParams_s(3:end);0;0];
   Y_ave=GlobalParams_s(:,2)+feckParams_s;
   Foods(:,2)=Y_ave(2*(r-node/2-runtime)-1);
   ObjVal=feval(objfun,Foods,GlobalParams_s(2*(r-node/2-runtime)-1,:),GlobalParams_s(2*(r-node/2-runtime)+1,:));            
    end    
end
%====================================================================================%

Fitness=calculateFitness(ObjVal);
%reset trial counters
trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal(find(ObjVal>0))));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
while (iter <= maxCycle),
% Employed bees phase - each employed bee visits a food source and searches its neighborhood
%%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        Param2Change=fix(rand*D)+1;
        neighbour=fix(rand*FoodNumber)+1;      
            while(neighbour==i) % Keep generating a random neighbor until it's different from i
                neighbour=fix(rand*FoodNumber)+1;
            end;
        
       sol=Foods(i,:);
       % Generate new solution using the current and neighbor solution
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2; % (rand-0.5)*2 generates a random number in [-1:1]
     
        ind=find(sol<(lb)); % If out of boundary, set to boundary value
        sol(ind)=lb(ind);
        ind=find(sol>(ub));
        sol(ind)=ub(ind);
        
        %evaluate new solution
        if r<=node/2
            if rem(r,2)~=0
          ObjValSol=feval(objfun,sol,GlobalParams_s(r,:),GlobalParams_s(node+4-r,:)); 
            else
          ObjValSol=feval(objfun,sol,GlobalParams_s(r+1,:),GlobalParams_s(node+5-r,:)); 
            end
      else
        if r>node/2&&r<=runtime
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2)-1,:),GlobalParams_s(2*(r-node/2)+1,:));   
        elseif r>runtime&&r<=runtime+node/2
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-runtime),:),GlobalParams_s(2*(r-runtime)+2,:));   
        else
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2-runtime)-1,:),GlobalParams_s(2*(r-node/2-runtime)+1,:));   
        end
       end
       
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant 
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
     end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob=(0.9.*Fitness./max(Fitness))+0.1; % ".*" denotes element-wise multiplication
%prob=Fitness./sum(Fitness);
%%%%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        Param2Change=fix(rand*D)+1; % Randomly select a parameter in range [1,D]
        neighbour=fix(rand*(FoodNumber))+1;
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol<(lb)); % If out of boundary, set to boundary value
        sol(ind)=lb(ind);
        ind=find(sol>(ub));
        sol(ind)=ub(ind);
        
        %evaluate new solution
        if r<=node/2
            if rem(r,2)~=0
          ObjValSol=feval(objfun,sol,GlobalParams_s(r,:),GlobalParams_s(node+4-r,:)); 
            else
          ObjValSol=feval(objfun,sol,GlobalParams_s(r+1,:),GlobalParams_s(node+5-r,:)); 
            end
      else
        if r>node/2&&r<=runtime
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2)-1,:),GlobalParams_s(2*(r-node/2)+1,:));   
        elseif r>runtime&&r<=runtime+node/2
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-runtime),:),GlobalParams_s(2*(r-runtime)+2,:));   
        else
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2-runtime)-1,:),GlobalParams_s(2*(r-node/2-runtime)+1,:));   
        end
       end
        
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
     if (i==(FoodNumber)+1)
         i=1;
     end;   
end; 


%/*The best food source is memorized*/
         ind=find(ObjVal==min(ObjVal(find(ObjVal>0))));
         ind=ind(end);
         if (ObjVal(ind)<GlobalMin)
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:); % Save parameters with smallest objective value
         end;
         
         
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used to prevent algorithm from getting stuck in local optima
ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    
    sol=(ub-lb).*rand(1,D)+lb; % Randomly generate new food source
    
     if r<=node/2
            if rem(r,2)~=0
          ObjValSol=feval(objfun,sol,GlobalParams_s(r,:),GlobalParams_s(node+4-r,:)); 
            else
          ObjValSol=feval(objfun,sol,GlobalParams_s(r+1,:),GlobalParams_s(node+5-r,:)); 
            end
      else
        if r>node/2&&r<=runtime
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2)-1,:),GlobalParams_s(2*(r-node/2)+1,:));   
        elseif r>runtime&&r<=runtime+node/2
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-runtime),:),GlobalParams_s(2*(r-runtime)+2,:));   
        else
        ObjValSol=feval(objfun,sol,GlobalParams_s(2*(r-node/2-runtime)-1,:),GlobalParams_s(2*(r-node/2-runtime)+1,:));   
        end
       end
     
    FitnessSol=calculateFitness(ObjValSol); % Calculate fitness
    Foods(ind,:)=sol; % Replace abandoned food source with newly generated food source
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;

fprintf('Iteration=%d GlobalMin=%g\n',iter,GlobalMin);
iter=iter+1;
end ; % End of ABC

GlobalMins(r)=GlobalMin;

if r<=node/2
    if rem(r,2)~=0
    GlobalParams_s(r+2,:)=GlobalParams;
    else
    GlobalParams_s(node+3-r,:)=GlobalParams;        
    end
else
  if r>node/2&&r<=runtime
  GlobalParams_s(2*(r-node/2),:)=GlobalParams;
  elseif r>runtime&&r<=runtime+node/2
  GlobalParams_s(2*(r-runtime)+1,:)=GlobalParams;  
  else
  GlobalParams_s(2*(r-node/2-runtime),:)=GlobalParams;    
  end
end

end; %end of runs

%===================================== Path Visualization ========================================%
% Choose a color for this UAV's path based on the obj parameter
colors = {'r', 'g', 'b', 'c', 'm', 'y'};
colorIdx = mod(obj-1, length(colors)) + 1;
uavColor = colors{colorIdx};

% Mark start and end points with distinctive markers
scatter3(setstart(1), setstart(2), setstart(3), 100, 'filled', 'MarkerFaceColor', uavColor, 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter3(setfinal(1), setfinal(2), setfinal(3), 100, 'filled', 'MarkerFaceColor', uavColor, 'MarkerEdgeColor', 'k', 'Marker', 's');

% Add labels for each UAV
text(setstart(1)+10, setstart(2)+10, setstart(3)+10, ['UAV #', num2str(obj)], 'Color', uavColor, 'FontWeight', 'bold');

% Plot the complete path
plot3(GlobalParams_s(:,1), GlobalParams_s(:,2), GlobalParams_s(:,3), 'Color', uavColor, 'LineWidth', 2);
scatter3(GlobalParams_s(:,1), GlobalParams_s(:,2), GlobalParams_s(:,3), 20, uavColor, 'filled');

% Add UAV height visualizer (vertical lines connecting path to ground)
try
    for i = 1:size(GlobalParams_s, 1)
        % Get the ground elevation at this point
        ground_z = interp2(X, Y, Z, GlobalParams_s(i,1), GlobalParams_s(i,2), 'cubic');
        if ~isnan(ground_z)
            % Draw vertical line from ground to flight path
            line([GlobalParams_s(i,1) GlobalParams_s(i,1)], ...
                 [GlobalParams_s(i,2) GlobalParams_s(i,2)], ...
                 [ground_z GlobalParams_s(i,3)], ...
                 'Color', uavColor, 'LineStyle', ':', 'LineWidth', 0.5);
        end
    end
catch
    warning('Could not draw height reference lines. Ground elevation data may be missing.');
end

% Save the calculated path data
save(['UAV_', num2str(obj), '_path.mat'], 'GlobalParams_s');

% Try to apply path smoothing visualization if the CBI function exists
try
    CBI; % Path smoothing visualization
catch
    warning('CBI function not found. Path smoothing visualization skipped.');
end

% Calculate path metrics
path_length = 0;
for i = 1:(size(GlobalParams_s,1)-1)
    segment_length = sqrt(sum((GlobalParams_s(i+1,:) - GlobalParams_s(i,:)).^2));
    path_length = path_length + segment_length;
end

% Calculate direct distance
direct_dist = sqrt(sum((setfinal - setstart).^2));

% Display path metrics in command window
fprintf('UAV #%d Path Metrics:\n', obj);
fprintf('  Direct distance: %.2f\n', direct_dist);
fprintf('  Path length: %.2f\n', path_length);
fprintf('  Path/Direct ratio: %.2f\n', path_length/direct_dist);

% Return the calculated path
POS = GlobalParams_s;

% Keep hold on so multiple UAVs can be visualized
hold on;