
%% Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing


close all
clear
wannaplot=1;
nFiles = 1998;%3768;

densityBL = [0.01,0.1];
densityAP = [9 12];%(1:1:10)/10^4;
omegaVal = pi/3;
connectivity = 1:4;
Arch = {'3GPP','FIBR'};

mu=2;

nMC = length(connectivity);
nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);
nArch = length(Arch);

tempInd = 0;
% tempInd2=0;
num0BS = zeros(1,nBL*nAP*nO);
num0BS_debug = zeros(1,nBL*nAP*nO);
Directory = {'3GPP_9/','3GPP_12/','FIBR_9/','FIBR_12/'};

prob_block = zeros(nArch,nAP,nMC,nBL);

prob_RLF = zeros(nArch,nAP,nMC,nBL);

block_dur = zeros(nArch,nAP,nMC,nBL);

thr = zeros(nArch,nAP,nMC,nBL);


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


for dir=1:4
    for indMC=1:nMC
        for indBL=1:nBL
            nonEmpty = 0;
            block_dur_greater_0 = 0;
            for aID=1:1998
                if (exist(strcat(Directory{dir},'output',int2str(connectivity(indMC)),'_',int2str(numBSs(dir)),'_',int2str(indBL),'_',int2str(aID),'.csv'))==0)
                    continue;
                else
                    %nonEmpty = nonEmpty + 1;
                    dataa = csvread(strcat(Directory{dir},'output',int2str(connectivity(indMC)),'_',int2str(numBSs(dir)),'_',int2str(indBL),'_',int2str(aID),'.csv'));
                    % Now dataa is a vector, with block_prob the 1st
                    % element, RLF_prob the 2nd, block_duration the 3rd,
                    % and the rest are the block times
                    %prob_block(archi(dir),iBS(dir),indMC,indBL) = prob_block(archi(dir),iBS(dir),indMC,indBL) + dataa(1);
                    %prob_RLF(archi(dir),iBS(dir),indMC,indBL) = prob_RLF(archi(dir),iBS(dir),indMC,indBL) + dataa(2);
                    if ~isnan(dataa(3)) % duration is greater than 0
                        block_dur_greater_0 = block_dur_greater_0 + 1;
                        prob_block(archi(dir),iBS(dir),indMC,indBL) = prob_block(archi(dir),iBS(dir),indMC,indBL) + dataa(1);
                        prob_RLF(archi(dir),iBS(dir),indMC,indBL) = prob_RLF(archi(dir),iBS(dir),indMC,indBL) + dataa(2);
                        block_dur(archi(dir),iBS(dir),indMC,indBL) = block_dur(archi(dir),iBS(dir),indMC,indBL) + dataa(3);                            
                    end % else just add nothing 
                end 
            end 
            prob_block(archi(dir),iBS(dir),indMC,indBL) = prob_block(archi(dir),iBS(dir),indMC,indBL)/block_dur_greater_0;
            prob_RLF(archi(dir),iBS(dir),indMC,indBL) = prob_RLF(archi(dir),iBS(dir),indMC,indBL)/block_dur_greater_0;
            block_dur(archi(dir),iBS(dir),indMC,indBL) = block_dur(archi(dir),iBS(dir),indMC,indBL)/block_dur_greater_0;
        end 
    end 
end

%% Process the above here
load('aver_Prob_block.mat');
load('RLF_prob.mat');
%load('block_dur.mat');
%load('prob_block.mat');
%load('prob_RLF.mat');
vector = zeros(nArch,2*nMC*2+1*nMC+2); % theoretical and simulation RLF and block probabilities for each MC degree, min RLF probability, min block probability, average blockaged urations for each MC degree 

for indBL=1:nBL
    for indBS=1:nAP
        vector = zeros(nArch,2*nMC*2+1*nMC+2);
        for indArch=1:nArch
            vector(indArch,1:4) = aver_Prob_block(indBL,1:4);
            vector(indArch,5:8) = RLF_prob(indBL,1:4);

            vector(indArch,21) = aver_Prob_block(indBL,4+indBS);
            vector(indArch,22) = RLF_prob(indBL,4+indBS);
            for indMC=1:nMC
                vector(indArch,8+indMC) = block_dur(indArch,indBS,indMC,indBL);
                vector(indArch,12+indMC) = prob_block(indArch,indBS,indMC,indBL);
                vector(indArch,16+indMC) = prob_RLF(indArch,indBS,indMC,indBL);
            end
        end 
        dlmwrite(strcat('plot','_blockerID',num2str(indBL),'_numBS',num2str(densityAP(indBS)),'.csv'),vector,'delimiter', ',','precision', 7)
    end 
end 


%% Throughput 

close all
clear
wannaplot=1;
nFiles = 1998;%3768;

densityBL = [0.01,0.1];
densityAP = [9 12];%(1:1:10)/10^4;
omegaVal = pi/3;
connectivity = 1:4;
Arch = {'3GPP','FIBR'};

mu=2;

nMC = length(connectivity);
nBL = length(densityBL);
nAP = length(densityAP);
nO = length(omegaVal);
nArch = length(Arch);

tempInd = 0;
% tempInd2=0;
num0BS = zeros(1,nBL*nAP*nO);
num0BS_debug = zeros(1,nBL*nAP*nO);
Directory = {'3GPP_9/','3GPP_12/','FIBR_9/','FIBR_12/'};

prob_block = zeros(nArch,nAP,nMC,nBL);

prob_RLF = zeros(nArch,nAP,nMC,nBL);

block_dur = zeros(nArch,nAP,nMC,nBL);

thr = zeros(nArch,nAP,nMC,nBL);


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
simTime = 4*60*60;


load('BS_distances_12.mat');
load('BS_distances_9.mat');
b = cell(4,1);
b{1} = BS_distances;
b{2} = BS_distances_12;
b{3} = BS_distances;
b{4} = BS_distances_12;
aver_throughput = zeros(nArch,nAP,nMC,nBL);

for dir=1:4
    for indMC=1:nMC
        for indBL=1:nBL
            nonEmpty = 0;
            through = 0;
            exists = 0;
            for aID=1:1998
                if (exist(strcat(Directory{dir},'output',int2str(connectivity(indMC)),'_',int2str(numBSs(dir)),'_',int2str(indBL),'_',int2str(aID),'.csv'))==0)
                    continue;
                else
                    exists = exists + 1;
                    dataa = csvread(strcat(Directory{dir},'output',int2str(connectivity(indMC)),'_',int2str(numBSs(dir)),'_',int2str(indBL),'_',int2str(aID),'.csv'));
                    if ~isnan(dataa(3)) % duration is greater than 0
                        nonEmpty = nonEmpty + 1;
                        
                        numBlock = (length(dataa) - 3)/2;
                        block_dur_inst = [dataa(4:4+numBlock-1)'];
                        block_time = dataa(4+numBlock:4+2*numBlock-1)';
                        
                        block_finish = [0 block_dur_inst + block_time];
                        on_Periods = [block_dur_inst - block_finish(1:end-1) simTime-block_finish(end)];
                        
                        through = through + distBasedThroughput((on_Periods),b{dir});
                    else % else just add nothing 
                        through = through + distBasedThroughput(simTime,min(b{dir}));
                    end 
                end 
            end
            % now we have found the average throughput for one aID
            aver_throughput(archi(dir),iBS(dir),indMC,indBL) = through/exists;
        end 
    end 
end

BS9_BL001 = [aver_throughput(1,1,1,1) aver_throughput(1,1,2,1) aver_throughput(1,1,3,1) aver_throughput(1,1,4,1); aver_throughput(2,1,1,1) aver_throughput(2,1,2,1) aver_throughput(2,1,3,1) aver_throughput(2,1,4,1)]
BS12_BL001 = [aver_throughput(1,2,1,1) aver_throughput(1,2,2,1) aver_throughput(1,2,3,1) aver_throughput(1,2,4,1); aver_throughput(2,2,1,1) aver_throughput(2,2,2,1) aver_throughput(2,2,3,1) aver_throughput(2,2,4,1)]
BS12_BL01 = [aver_throughput(1,2,1,2) aver_throughput(1,2,2,2) aver_throughput(1,2,3,2) aver_throughput(1,2,4,2); aver_throughput(2,2,1,2) aver_throughput(2,2,2,2) aver_throughput(2,2,3,1) aver_throughput(2,2,4,2)]
BS9_BL01 = [aver_throughput(1,1,1,2) aver_throughput(1,1,2,2) aver_throughput(1,1,3,2) aver_throughput(1,1,4,2); aver_throughput(2,1,1,2) aver_throughput(2,1,2,2) aver_throughput(2,1,3,1) aver_throughput(2,1,4,2)]


dlmwrite('thr_BS9_BL001.csv',BS9_BL001,'delimiter', ',', 'precision', 7)
dlmwrite('thr_BS12_BL001.csv',BS12_BL001,'delimiter', ',', 'precision', 7)
dlmwrite('thr_BS9_BL01.csv',BS9_BL01,'delimiter', ',', 'precision', 7)
dlmwrite('thr_BS12_BL01.csv',BS12_BL01,'delimiter', ',', 'precision', 7)

