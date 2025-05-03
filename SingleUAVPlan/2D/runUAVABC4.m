clc
close all
clear all
%===================================== Parameter Settings ==========================================%
global boundary setstart setfinal node danger_xi danger_yi danger_ri weight setstart_2 setfinal_2 ;
%node=12;
danger_xi=[250 100 300 50 400 170 156];% Obstacle center x-coordinates
danger_yi=[150 100 350 300 100 300 207];
%danger_ri=[80 50 80];% Obstacle radii
danger_ri=[60 65 70 50 60 30 20];
weight=[1 0.9];% Path length and collision risk weights
boundary=[500 0];% Environment boundaries
setfinal=[325 150];% End point
setstart=[0 230];% Start point
L_ALL=sqrt(((setfinal-setstart).^2)*ones(2,1));
L_FEN=30;% Node spacing distance
node=floor(L_ALL/L_FEN);
if mod(node,2)==1
    node=node+1;
end;
%===================================== Obstacle Visualization =======================================%
circle(danger_xi,danger_yi,danger_ri);
%========================================================================================%
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
            while(neighbour==i)% Keep generating random numbers until different from i
                neighbour=fix(rand*FoodNumber)+1;
            end;
        
       sol=Foods(i,:);
       % Generate new solution using the current and neighbor solution
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;% (rand-0.5)*2 generates a random number in [-1:1]
      % sol(1)=Foods(i,1);
     
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
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
     end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob=(0.9.*Fitness./max(Fitness))+0.1;% "." indicates element-wise operation
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
            end;
        
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
    
    sol=(ub-lb).*rand(1,D)+lb;% Randomly generate new food source
  
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
     
        
    FitnessSol=calculateFitness(ObjValSol);% Fitness calculation
    Foods(ind,:)=sol;% Food sources visited over limit times are replaced by newly generated food sources
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
plot([setstart(1) setfinal(1)],[setstart(2) setfinal(2)],'*');
%line([setfinal(1) GlobalParams_s(node+2,1) ],[setfinal(2) GlobalParams_s(node+2,2)],'Color','k','LineWidth',1);
%line([setstart(1) GlobalParams_s(1,1) ],[setstart(2) GlobalParams_s(1,2)],'Color','k','LineWidth',1);
%plot(GlobalParams_s(:,1),GlobalParams_s(:,2),'*','Color','r');
%line(GlobalParams_s(:,1),GlobalParams_s(:,2),'Color','k','LineWidth',1); 
%====================================================================================%
B;% B-spline smoothing
hold on;
save all
%================================ Path Distance Calculation =========================================%
C=GlobalParams_s;
A=0;
for i=1:(size(C,1)-1)
A=A+sqrt(((C(i+1,:)-C(i,:)).^2)*ones(2,1));
end
min=sqrt(((setstart-setfinal).^2)*ones(2,1))% Direct distance from start to end
routes_distance=A;
routes_distance% Path planning total distance