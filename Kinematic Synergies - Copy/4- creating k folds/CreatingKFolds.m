% This mfile reads the fullset data for a set of subjects and creates k folds of
% trainging and validation data for each subject.


% 20151209 Written by Navid Shirzad



load('FullSet_SS3.mat');
RightFullSet = [FullSet(2:end,2:6) FullSet(2:end,8:12)]; %10DOFs
LeftFullSet = [FullSet(2:end,2:4) FullSet(2:end,13:14) FullSet(2:end,16:20)]; %2:end because I don't want to include t=0s

% divide the data into 5s epochs
numEpochs = 315/5;

for k=1:20 % k is number of folds
    RightTrain = zeros(120,10);
    LeftTrain = zeros(120,10);
    RightVal = zeros(30,10);
    LeftVal = zeros(30,10);

    for temp=1:numEpochs
       %in each epoch there are 150 data points , 20% of the data will be
       %assigned to validation (30 data points) and the rest for training.
       %Choose the start of the 30 data points randomly
       StartVal = randi([2,149]);
       if StartVal<121
           RightEpochVal = RightFullSet((temp-1)*150+StartVal:(temp-1)*150+StartVal+29,:);
           LeftEpochVal = LeftFullSet((temp-1)*150+StartVal:(temp-1)*150+StartVal+29,:);
           RightEpochTr = [RightFullSet((temp-1)*150+1:(temp-1)*150+StartVal-1,:); RightFullSet((temp-1)*150+StartVal+30:temp*150,:)];
           LeftEpochTr = [LeftFullSet((temp-1)*150+1:(temp-1)*150+StartVal-1,:); LeftFullSet((temp-1)*150+StartVal+30:temp*150,:)];
       else
           RightEpochTr = RightFullSet((temp-1)*150+StartVal-120:(temp-1)*150+StartVal-1,:);
           LeftEpochTr = LeftFullSet((temp-1)*150+StartVal-120:(temp-1)*150+StartVal-1,:);
           RightEpochVal = [RightFullSet((temp-1)*150+1:(temp-1)*150+StartVal-121,:); RightFullSet((temp-1)*150+StartVal:temp*150,:)];
           LeftEpochVal = [LeftFullSet((temp-1)*150+1:(temp-1)*150+StartVal-121,:); LeftFullSet((temp-1)*150+StartVal:temp*150,:)];
       end

       if temp == 1 
           RightTrain = RightEpochTr;
           LeftTrain = LeftEpochTr;
           RightVal = RightEpochVal;
           LeftVal = LeftEpochVal;
       else          
           RightTrain = [RightTrain; RightEpochTr];
           LeftTrain = [LeftTrain; LeftEpochTr];
           RightVal = [RightVal; RightEpochVal];
           LeftVal = [LeftVal; LeftEpochVal];
       end

    end
    
    save(strcat('SS3_Right_Tr_Fold_', num2str(k), '.mat' ), 'RightTrain')
    save(strcat('SS3_Left_Tr_Fold_', num2str(k), '.mat' ), 'LeftTrain')
    save(strcat('SS3_Right_Val_Fold_', num2str(k), '.mat' ), 'RightVal')
    save(strcat('SS3_Left_Val_Fold_', num2str(k), '.mat' ), 'LeftVal')
    
end