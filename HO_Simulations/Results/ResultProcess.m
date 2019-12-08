
%% Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing


close all
clear
wannaplot=1;
nFiles = 4000;


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
densityAP = poissrnd(densityBS*pi*R^2);
connectivity = [1 2 3 4 ];

omegaVal = pi/3;

mu=2;

nMC = length(connectivity);
nBL = length(densityBL);
nAP = length(densityAP);
nDisc = length(discovery);
nPrep = length(preparation);



tempInd = 0;
% tempInd2=0;
num0BS = zeros(1,nBL*nAP*nO);
num0BS_debug = zeros(1,nBL*nAP*nO);
Directory = {'3GPP_9/','3GPP_12/','FIBR_9/','FIBR_12/'};

prob_block = zeros(nAP,nMC,nBL,nPrep,nDisc);

prob_RLF = zeros(nAP,nMC,nBL,nPrep,nDisc);

block_dur = zeros(nAP,nMC,nBL,nPrep,nDisc);

thr = zeros(nAP,nMC,nBL,nPrep,nDisc);


numBSs = zeros(1,length(Directory));
iBS(1) = 1; %  1 for 9, 2 for 12, useful for indexing
iBS(2) = 2;
iBS(3) = 1;
iBS(4) = 2;

numBSs(1) = 9; % now this is useful for the name 
numBSs(2) = 12;
numBSs(3) = 9;
numBSs(4) = 12;

% archi is 1 for 3GPP, 2 for FIBR
archi = zeros(1,length(Directory));
archi(1) = 1;
archi(2) = 1;
archi(3) = 2;
archi(4) = 2;



for aID=1:6000
    %nonEmpty = nonEmpty + 1;
    dataa = csvread(strcat('output',int2str(aID),'.csv'));
    for rowData=1:10 % we have 30 columns but 3 consecutive elements belong to the same ind
        for colData=1:32
            % First we find the corresponding indeces
            indAP = ceil(colData/8); % number of BS elements indeces are repeated every 8 values
            indMC = ceil(mod(colData,8)/2);
            if indMC == 0
                indMC = 4; 
            end
            indBD = mod(colData,2); %number of block densities elements indeces keep inter-changing
            if indBD == 0
                indBD = 2;
            end
            indPrep = mod(rowData,2);
            if indPrep == 0
                indPrep = 2;
            end
            indDisc = ceil(rowData/2);
            % Now, we update probabilities and blockage duration
            prob_block(indAP,indMC,indBL,indPrep,indDisc) = prob_block(indAP,indMC,indBL,indPrep,indDisc) + data(rowData,colData);
            prob_RLF(indAP,indMC,indBL,indPrep,indDisc) = prob_RLF(indAP,indMC,indBL,indPrep,indDisc) + data(rowData+1,colData);
            block_dur(indAP,indMC,indBL,indPrep,indDisc) = block_dur(indAP,indMC,indBL,indPrep,indDisc) + data(rowData+2,colData);
        end
    end    
end 
prob_block = prob_block/6000;
prob_RLF = prob_RLF/6000;
block_dur = block_dur/6000;

%% Process the above here

% need to load sthe theoretical results here to compare with the simulation


%load('aver_Prob_block.mat');


for indBL=1:nBL
    for indBS=1:nAP
        % need to add stuff here based on what we want to show
        dlmwrite(strcat('plot','_blockerID',num2str(indBL),'_numBS',num2str(densityAP(indBS)),'.csv'),vector,'delimiter', ',','precision', 7)
    end 
end 


