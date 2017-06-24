% This mfile reads the fullset data for a set of subjects and creates k folds of
% trainging and validation data for each subject. 
% try running CreatingKFolds_BatchProcess([1 5], 20) to create 20 folds for
% each of the two subjects, i.e., #1 and #5.


% 20151209 Written by Navid Shirzad



function CreatingKFolds_BatchProcess(SubjectIDs, NumberFolds)

    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end

        % Read in the data 
        load(strcat('FullSet_', 'Y', SubjID, '.mat'))
        RightFullSet = [FullSet(:,2:6) FullSet(:,8:12)]; %10DOFs
        LeftFullSet = [FullSet(:,2:4) FullSet(:,13:14) FullSet(:,16:20)]; %t=0s is not included in FullSet (first row is t=1/30s)

        % divide the data into 5s epochs
        numEpochs = FullSet(end,1)/5;

        for k=1:NumberFolds % k is number of folds
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

            save(strcat('Y', SubjID, '_Right_Tr_Fold_', num2str(k), '.mat' ), 'RightTrain')
            save(strcat('Y', SubjID, '_Left_Tr_Fold_', num2str(k), '.mat' ), 'LeftTrain')
            save(strcat('Y', SubjID, '_Right_Val_Fold_', num2str(k), '.mat' ), 'RightVal')
            save(strcat('Y', SubjID, '_Left_Val_Fold_', num2str(k), '.mat' ), 'LeftVal')

        end
    
    end
end