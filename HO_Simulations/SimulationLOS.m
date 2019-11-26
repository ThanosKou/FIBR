% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Description:
% First we get the blocker mobility using Generate_Mobility.m function.
% Then for different BS Densities, blocker densities and self-blockage
% angle, we call BlockageSimFn.m function to get key blockage metrics like
% blockage duration, frequency, and blockage. We should run this code for
% many iterations prefebly on high performance computing machine.

close all;
clear;

%----Play-with-values---------------------------------------
aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
rng('shuffle');

% considerLOS=0;
% considerNLOS=1;
wannaplot=0; %Don't plot if running for many loops (else too many plots).
V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
simTime = 4*60*60; %sec Total Simulation time
% Note!!! simTime must be >100s else the code won't work :)
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius

discovery = [1 5 20]*10^(-3);
preparation = [10 20]*1^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4 ];

nTorig = densityBS*pi*R^2;
omega = pi/3;

s_input = cell(1,2);
s_mobility = cell(1,2);

for indB=1:length(densityBL)
s_input{indB} = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[V V],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
    'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',4*R^2*densityBL(indB));

% Generate_Mobility function is Copyright (c) 2011, Mathieu Boutin
s_mobility{indB} = Generate_Mobility(s_input{indB});
end
%finaldata = zeros(3,length(connectivity),length(densityBL),length(discovery),length(preparation),length(densityBS));

for indBS = 1:length(densityBS)
    nT = poissrnd(densityBS(indBS)*pi*R^2);
    rT = R*sqrt(rand(nT,1)); %location of APs (distance from origin)
    alphaT = 2*pi*rand(nT,1);%location of APs (angle from x-axis)
    for indT = 1:length(connectivity)
        currConnec = connectivity(indT);
        for indB = 1:length(densityBL) %for all blockers
            rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
            nB = 4*R^2*rhoB;%=4000; %number of blokers

            BS_input = struct('WANNAPLOT',wannaplot,...
                'DEGREE_CONNECTIVITY', currConnec,...
                'RADIUS_AROUND_UE',R,...
                'SIMULATION_TIME',simTime,...
                'TIME_STEP',tstep,...
                'MU',mu,...
                'FRACTION',frac,...
                'SELF_BL_ANGLE_OMEGA',omega,...
                'Original_NUM_AP',nT,...
                'LOC_AP_DISTANCE', rT,...
                'LOC_AP_ANGLE',alphaT,...
                'NUM_BL',nB);

                %BlockageSimFn function is written by Ish Jain
                [BS_state,output] = BlockageSimFn(s_mobility{indB},BS_input);
                discovery = [1 5 20];
                preparation = [10 20];
                indCell = 0; 
                for disc=1:length(discovery)
                    for prep=1:length(preparation)
                        indCell = indCell + 1;
                        dlmwrite(strcat('output',num2str(nT),'_',num2str(currConnec),'_',num2str(indB),'_',num2str(discovery(disc)),'_',num2str(preparation(prep)),'_',num2str(aID),'.csv'),output{indCell},'delimiter',',','precision',7)
                        dlmwrite(strcat('BS_state',num2str(nT),'_',num2str(currConnec),'_',num2str(indB),'_',num2str(discovery(disc)),'_',num2str(preparation(prep)),'_',num2str(aID),'.csv'),BS_state{indCell},'delimiter',',','precision',7)
                    end 
                end 
                %finaldata(:,indBS,indT,indB) = output(1:3);
            %         output is [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
        end
   end
end 
%Use the code processData9.m to analyze and plot the results
%csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)
