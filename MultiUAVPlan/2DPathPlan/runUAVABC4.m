function POS=runUAVABC4(obj)
global boundary setstart setfinal node setstart_2 setfinal_2 ;
maxCycle=10;% Maximum algorithm iterations
%/* Problem specific variables*/

D=2; 
NP=40; 
FoodNumber=NP/2; 
limit=NP*D;
ub=ones(1,D)*(boundary(1)-10); 
lb=ones(1,D)*(boundary(2)+10);
setfinal_2=setfinal;
setstart_2=setstart;
runtime=node+1;
objfun='UAV2';

GlobalMins=zeros(1,runtime);
% Fix for matrix dimension error - ensure setstart_2 and setfinal_2 are row vectors
if size(setstart_2,1) > 1
    setstart_2 = setstart_2';
end
if size(setfinal_2,1) > 1
    setfinal_2 = setfinal_2';
end

% Check if vectors have the correct number of columns
if size(setstart_2,2) ~= 2
    setstart_2 = setstart_2(1:2);
end
if size(setfinal_2,2) ~= 2
    setfinal_2 = setfinal_2(1:2);
end

GlobalParams_s=[setstart_2;zeros(node+1,2);setfinal_2];

for r=1:2*runtime
Range = repmat((ub-lb),[FoodNumber 1]);
Lower = repmat(lb, [FoodNumber 1]);
Foods = rand(FoodNumber,D) .* Range + Lower;

%====================================================================================%
if r<=node/2% Step1: Plan odd position nodes in pairs, distributing by X coordinates
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
    if r>node/2&&r<=runtime% Step2: Plan even position nodes individually, distributing by Y coordinates based on nodes from Step1
   feckParams_s=[0;0;GlobalParams_s(1:(end-2),2)];
   feckParams_s=(GlobalParams_s(:,2)-feckParams_s)/2;
   feckParams_s=[feckParams_s(3:end);0;0];
   Y_ave=GlobalParams_s(:,2)+feckParams_s;
   Foods(:,2)=Y_ave(2*(r-node/2)-1);
   ObjVal=feval(objfun,Foods,GlobalParams_s(2*(r-node/2)-1,:),GlobalParams_s(2*(r-node/2)+1,:));   
    elseif r>runtime&&r<=runtime+node/2% Step3: Refine odd position nodes individually, distributing by X coordinates based on nodes from Step2
   feckParams_s=[0;GlobalParams_s(4:end,1);0;0];     
   feckParams_s=(feckParams_s-GlobalParams_s(:,1))/2;  
   X_ave=GlobalParams_s(:,1)+feckParams_s;
   Foods(:,1)=X_ave(2*(r-runtime));
   ObjVal=feval(objfun,Foods,GlobalParams_s(2*(r-runtime),:),GlobalParams_s(2*(r-runtime)+2,:));   
    else% Step4: Refine even position nodes individually, distributing by Y coordinates based on nodes from Step3
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
            while(neighbour==i)
                neighbour=fix(rand*FoodNumber)+1;
            end;
        
       sol=Foods(i,:);
       % Generate new solution using the current and neighbor solution
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;% (rand-0.5)*2 generates a random number in [-1:1]
     
        ind=find(sol<(lb));% If out of boundary, set to boundary value
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
            trial(i)=trial(i)+1;
       end;
     end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob=(0.9.*Fitness./max(Fitness))+0.1;
%prob=Fitness./sum(Fitness);
%%%%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        Param2Change=fix(rand*D)+1;% Randomly select a parameter in range [1,D]
        neighbour=fix(rand*(FoodNumber))+1;
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end
        
       sol=Foods(i,:);
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol<(lb));% If out of boundary, set to boundary value
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
         GlobalParams=Foods(ind,:);% Save parameters with smallest objective value
         end;
         
         
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used to prevent algorithm from getting stuck in local optima
ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    
    sol=(ub-lb).*rand(1,D)+lb;
  
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
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;

fprintf('Iteration=%d GlobalMin=%g\n',iter,GlobalMin);
iter=iter+1;
end ;% End of ABC

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
% Mark start and end points
plot(setstart(1), setstart(2), 'b*', 'MarkerSize', 8);
plot(setfinal(1), setfinal(2), 'g*', 'MarkerSize', 8);

% Add UAV number label if provided
if nargin > 0 && ~isempty(obj)
    text(setstart(1)+10, setstart(2)+10, ['UAV #', num2str(obj)]);
end

% Plot path
plot(GlobalParams_s(:,1), GlobalParams_s(:,2), 'k-', 'LineWidth', 1.5);
plot(GlobalParams_s(:,1), GlobalParams_s(:,2), 'ro', 'MarkerSize', 4);

% Apply B-spline smoothing if available
try
    % Try running B-spline smoothing
    P = GlobalParams_s'; 
    n = size(P,2)-1; 
    k = 2; % B-spline degree
    NodeVector = U_quasi_uniform(n, k);
    DrawSpline(n, k, P, NodeVector);
    hold on;
catch
    % If B-spline fails, no need to do anything as we already plotted the path
end

hold on;
save all

%================================ Path Distance Calculation =========================================%
% Calculate path length
C = GlobalParams_s;
path_length = 0;
for i = 1:(size(C,1)-1)
    % Calculate segment length using Euclidean distance
    segment_length = sqrt(sum((C(i+1,:) - C(i,:)).^2));
    path_length = path_length + segment_length;
end

% Calculate direct distance from start to end
% Make sure setstart and setfinal are converted to same format for calculation
if size(setstart, 1) > 1
    start_point = setstart';
else
    start_point = setstart;
end

if size(setfinal, 1) > 1
    end_point = setfinal';
else
    end_point = setfinal;
end

direct_dist = sqrt(sum((end_point - start_point).^2));

% Display distance information
fprintf('Direct distance: %f\n', direct_dist);
fprintf('Path planning total distance: %f\n', path_length);

% Return final path
POS = GlobalParams_s;