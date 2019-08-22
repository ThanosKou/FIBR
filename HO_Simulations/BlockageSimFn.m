function output = BlockageSimFn(s_mobility,BS_input)
% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Input:
%   s_mobility: contains parameters corresponding to blockers mobility
%   BS_input: contains infromation related to BS-UE topology and simulation
%   parameters. See SimulationLOS.m for the usage.
% Description:
% We use random way point mobility model for the blockers. The UE is at the
% origin and the BSs are located in a circle around the UE. A random
% blocker will block the BS-UE LOS path is it crosses the line joining
% between BS and the UE. We use find_blockage_distance.m function to fine
% intersection of two such lines. Finally, we repeat the process for all
% the blockers, for all the BSs and for the whole simulation duration. We
% collect this data as sequences (binary_seq) of 0 and 1 (0=blockage, 1=not blocked) for
% the entire duration for each BS-UE pair. Bitwise and of these sequences
% for all the APs gives a sequence (allBl) indicating simultaneous blockage of all BSs.
%
% Output: output file contains:
%   avgFreq: Average frequency of simultaneous blockage of all BS. 
%   avgDur: Average duration of simultaneous blockage of all BS.
%   probAllBl: probability of simultaneous blockage of all BS.
%   nTorig: original number of BS (can be blocked by self-blockage)
%   nT: number of BS not blocked by self-blockage



%----Play-with-values-here--------------------------------------
wannaplot = BS_input.WANNAPLOT; %1;
nB = BS_input.NUM_BL; %number of blokers
nTorig = BS_input.Original_NUM_AP; %Original APs without considering self blockage
rT =BS_input.LOC_AP_DISTANCE; %location of APs
alphaTorig = BS_input.LOC_AP_ANGLE;%location of APs

frac = BS_input.FRACTION;
omega = BS_input.SELF_BL_ANGLE_OMEGA;  

%%Implementing self-blockage
tempInd =  find(alphaTorig>=omega); %These BSs are not blocked by self-blockage
xT = rT(tempInd).*cos(alphaTorig(tempInd));%location of APs (distance)
yT = rT(tempInd).*sin(alphaTorig(tempInd));%location of APs (angle)
nT = length(tempInd); % number of BS not blocked by self-blockage
% nT=0
if(nT==0)
    output=[0,0,0];   
    return;
end % Dealing zero APs
   
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns
alphaT = alphaTorig(tempInd); %angle from x-axis for BS not blocked by self-bl
simTime = BS_input.SIMULATION_TIME; %sec Total Simulation time
tstep = BS_input.TIME_STEP; %(sec) time step
mu = BS_input.MU; %Expected bloc dur =1/mu
conDegree = BS_input.DEGREE_CONNECTIVITY;

RLF_timer = 0.03;
RLF_recovery = 0.04;

dataBS = cell(nT,1); 
%dataBS contain array of timestamps of blocker arrival for all BSs,
%independent of which blocker is blocking

for indB = 1:nB %for every blocker
    
    for iter =1:(length(s_mobility.VS_NODE(indB).V_POSITION_X)-1)
        
        % for every time blocker changes direction
        loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter)];
        loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
        start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
        velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
        for indT = 1:nT %for every BS around the UE (outside self-bl zone)
            %The find_blockage_distance() function is written by Ish Jain
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),alphaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                dataBS{indT} = [dataBS{indT}, timestampBl];
                
            end
            
        end
        
    end
end


for i=1:nT
    dataBS{i} = sort(dataBS{i});
end

totaltime = (simTime)/tstep;
binary_seq = zeros(nT,totaltime); %time sequence for every BS
allBl = ones(1,totaltime); %binary seq of all blocked
Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
if(wannaplot),figure; hold on; end

tranBSs=1:1:nT;
BSSET = randperm(nT,conDegree);
NONBSSET = setdiff(tranBSs,BSSET);
BLOCKEDBSSET = [];

for indT = 1:nT
    len =length(dataBS{indT});
    dataBS{indT}(2,:) =  exprnd(1/mu,1,len); % block duration
end

%% Thanos & Rajeev, FIBR work

dt = 0.1; % time after getting blocked

servBS = zeros(4,nT);

choose = @(samples) samples(randi(numel(samples)));

tt = ones(1,nT); %loop index for all BSs
timestamp = 0;
actions = [];
blockage_duration = [];

while timestamp < totaltime
    
    if ~isempty(actions)
        l = struct2cell(actions);
        [x,action_index] = min(cell2mat(l(1,:))); % find next action 
        if x < timestamp    
            g = 1;
            timestamp = totaltime;
            continue
        else
                timestamp = x;
        end 
        %timestamp = actions(action_index).timeinstance;  % timestamp of next action
     
        if strcmp(actions(action_index).fnc,'add')
            if ~isempty(NONBSSET) % there are more BSs in the nonBSset          
                newBS = choose(NONBSSET); % pick a random BS and add it to the set 
                % Need to check if this BS is currently blocked
                last_time_index = find(dataBS{newBS}(1,:)<=timestamp,1,'last');
                last_time_blocked = dataBS{newBS}(1,last_time_index);
                last_blockage_dur = dataBS{newBS}(2,last_time_index);
                if last_time_blocked + last_blockage_dur >= timestamp
                    % it means that this BS is blocked, need to try
                    % again in the next dt
                    BLOCKEDBSSET = [BLOCKEDBSSET newBS];
                    NONBSSET = setdiff(NONBSSET,newBS);
                    actions = [actions struct('timeinstance',{timestamp + dt},'BSindex',{1},'fnc',{'add'})]; % add a new BS to BSSET
                    %actions = [actions struct('timeinstance',{dataBS{newBS}(1,tt(newBS))},'BSindex',{newBS},'fnc',{'nextBlock'})];
                    actions = [actions struct('timeinstance',{last_time_blocked + last_blockage_dur},'BSindex',{newBS},'fnc',{'recover'})];  % add it again to NONBSSET when blockage ends
                else
                    if isempty(BSSET) 
%                         if (timestamp - blockage_duration(end)  < RLF_timer) % if it was empty and now a new BS is added, then it stops being blocked
%                             blockage_duration(end) = timestamp - blockage_duration(end);
%                         end
                       blockage_duration(end) = timestamp - blockage_duration(end); 
                    end 
                    BSSET = [BSSET newBS]; 
                    NONBSSET = setdiff(NONBSSET,newBS);
                    
                    % need to take care of the possibility that the new BS
                    % will be immediately blocked
                    if ~isempty(find(dataBS{newBS}(1,:)>=timestamp,1,'first'))
                        tt(newBS) = find(dataBS{newBS}(1,:)>=timestamp,1,'first');
                        actions = [actions struct('timeinstance',{dataBS{newBS}(1,tt(newBS))},'BSindex',{newBS},'fnc',{'nextBlock'})];
                    end 
                end
            end
            
            actions(action_index) = []; % remove the current add action
            
            if isempty(BSSET) && (timestamp - blockage_duration(end)  >= RLF_timer)  % an RLF is declared after 30ms and we have to wait for another 40ms for RLF recovery procedure
            %if timestamp - blockage_duration(end)  >= RLF_timer  % an RLF is declared after 30ms and we have to wait for another 40ms for RLF recovery procedure
                actions = []; % clean all the actions and choose new BSs after 40ms
                timestamp = timestamp + RLF_recovery;
                
                BSSET = randperm(nT,conDegree);
                NONBSSET = setdiff(tranBSs,BSSET);
                BLOCKEDBSSET = [];
                
                % Now that we randonly selected the BSs that we
                % start with, we should check if they will be up
                % after 40ms.
                BSafterRecovery = 0 ;
                BSSET_copy = BSSET;
                for indBS = 1:length(BSSET)
                    last_time_index = find(dataBS{BSSET(indBS)}(1,:)<=timestamp,1,'last');
                    last_time_blocked = dataBS{BSSET(indBS)}(1,last_time_index);
                    last_blockage_dur = dataBS{BSSET(indBS)}(2,last_time_index);
                    if last_time_blocked + last_blockage_dur >= timestamp
                        actions = [actions struct('timeinstance',{last_time_blocked + last_blockage_dur},'BSindex',{BSSET(indBS)},'fnc',{'recover'})];  % add it again to NONBSSET when blockage ends
                        actions = [actions struct('timeinstance',{timestamp + dt},'BSindex',{1},'fnc',{'add'})];
                        old_bs = BSSET(indBS);
                        BLOCKEDBSSET = [BLOCKEDBSSET old_bs]; 
                        BSSET_copy = setdiff(BSSET_copy,old_bs);
                    else
                        BSafterRecovery = BSafterRecovery + 1;  % we may not need this value, but if it is greater than 0, then we will have connectivity after RLF recovery 
                    end         
                end 
                BSSET = BSSET_copy;
                if BSafterRecovery > 0
                    blockage_duration(end) = timestamp - blockage_duration(end); % the blockage is ended!
                    actions = [actions struct('timeinstance',{timestamp},'BSindex',{100},'fnc',{'nextBlock'})]; % this is a trick to bypass next iteration actions 'add' and 'recover' 
                else
                    actions = [actions struct('timeinstance',{timestamp + dt + RLF_recovery},'BSindex',{1},'fnc',{'add'})]; % add a new BS to BSSET    
                end
            end 
            
            continue
        elseif strcmp(actions(action_index).fnc,'recover')
            recoveredBS = actions(action_index).BSindex;
            NONBSSET = [NONBSSET recoveredBS];
            BLOCKEDBSSET = setdiff(BLOCKEDBSSET,recoveredBS);
            actions(action_index) = [];
            continue
        else 
            actions(action_index) = [];
        end 
    end 
    
    % First, we get the blockage times for our serving BSs
    if isempty(BSSET)
        if isempty(actions) 
            actions = [actions struct('timeinstance',{timestamp + dt},'BSindex',{1},'fnc',{'add'})]; % add a new BS to BSSET
        else
            continue
        end 
    else 
        finishedBS = 0; % if finishedBS == nT, I am done
        for indBS = 1:nT
            if any(BSSET(:) == indBS) % if the current BS serves the UE
                %if timestamp > dataBS{BSSET(indBS)}(1,tt(indBS)) 
                if ~isempty(find(dataBS{indBS}(1,:)>=timestamp,1,'first'))
                    tt(indT) = find(dataBS{indBS}(1,:)>=timestamp,1,'first'); % find first blockage time index after the current timestamp th

                    servBS(1,indBS) = dataBS{indBS}(1,tt(indT)); % blockage times for serving BS
                    servBS(2,indBS) = dataBS{indBS}(2,tt(indT)); % blockage duration
                    servBS(3,indBS) = servBS(1,indBS) + dt; % UE cannot change BS until that time
                    servBS(4,indBS) = servBS(1,indBS) + servBS(2,indBS); % BS is again up at this time
                else
                    finishedBS = finishedBS + 1;
                %   timestamp = totaltime;
                %  continue
                end 
            %else
            %    servBS(:,indBS) = 0;
            end 
        end
        
        if finishedBS == length(BSSET)
            timestamp = totaltime;
            continue
        end 
        

        % Then, we find which BS is blocked first
        A = servBS(1,:);        
        if ~isempty(A(BSSET))
            [blockTime,blockedBS] = min(A(BSSET));
            old_bs = BSSET(blockedBS);
            BLOCKEDBSSET = [BLOCKEDBSSET old_bs]; 
            BSSET = setdiff(BSSET,old_bs); % remove the blocked BS from the list and then add it to blocked set
            
            
            if ~isempty(A(BSSET))
                [blockTimeNext,blockedBSNext] = min(A(BSSET)); % Prepare for the next blocked BS
                if ~isempty(actions)
                    ll = struct2cell(actions);
                    actions = actions(find(~strcmp(cellstr(ll(3,:)),'nextBlock')));
                end    
                actions = [actions struct('timeinstance',{blockTimeNext},'BSindex',{100},'fnc',{'nextBlock'})]; % next BS to be blocked, if the new ones do not get blocked before
            else
                blockage_duration = [blockage_duration timestamp];
            end
            actions = [actions struct('timeinstance',{servBS(3,old_bs)},'BSindex',{old_bs},'fnc',{'add'})]; % add a new BS to BSSET
            actions = [actions struct('timeinstance',{servBS(4,old_bs)},'BSindex',{old_bs},'fnc',{'recover'})];  % add it again to NONBSSET when blockage ends
        else  % that means that the UE just got disconnected from all BSs and I can start counting 
            %blockage_duration = [blockage_duration timestamp];
            
        end 

        %actions = [actions struct('timestamp',{blockTime},'BSindex',{blockedBS},'function',{'remove'})];

        % The next two actions need to be repeated every time a BS is blocked

    end 
end 

%%Evaluate frequency and average duration of blockage
%%frequency = (per sec)count transitions from 0 to 1 devided by simTime
%%duration: (sec) count total # of 1's multiplied by tstep/blockageCount
% 
% avgFreq = sum(diff(allBl)>0)/simTime;
% avgDur = sum(allBl)*tstep/sum(diff(allBl)>0);
% probAllBl = sum(allBl)*tstep/simTime;


avgDur = mean(blockage_duration(2:end));
probBl = sum(blockage_duration(2:end))/simTime;
RLFprob = sum(blockage_duration(blockage_duration(2:end)>0.03))/simTime;
%avgFreq = length(blockage_duration)/simTime;

%%Return now
output=[probBl,RLFprob,avgDur]';

end
