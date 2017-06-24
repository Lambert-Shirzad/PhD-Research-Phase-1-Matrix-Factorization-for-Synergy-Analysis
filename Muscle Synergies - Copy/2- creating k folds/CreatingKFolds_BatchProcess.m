% This mfile reads the filtered+RMS+Normalized EMG data stored in separate
% mat files for right and left side of thebody and creates k folds of
% trainging and validation data for each subject. 
% try running CreatingKFolds_BatchProcess([1 5], 20) to create 20 folds for
% each of the two subjects, i.e., #1 and #5.

%change lines 25,26,73:76 to go between Weighted and Lucky Pirate versions of
%the data.

% 20160322 Written by Navid Shirzad



function CreatingKFolds_BatchProcess(SubjectIDs, NumberFolds)

    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end

        % Read in the data 
        load(strcat('Processed_Weight_Subj_', SubjID, '_Right.mat'))
        load(strcat('Processed_Weight_Subj_', SubjID, '_Left.mat'))
        
        % divide the data into 5s epochs      
        numEpochs = floor(ProcessedRightSide(end,1)/5); %look at the time, how many 5s epochs are there?
        if numEpochs > 300/5
            numEpochs = 300/5; % originally we have about 305s of data, use 300s of the data %ProcessedRightSide(end,1)/5;
        end
        ProcessedRightSide = ProcessedRightSide(:,2:9); %8EMG Channels
        ProcessedLeftSide = ProcessedLeftSide(:,2:9); %ignore column 1 (time)

        for k=1:NumberFolds % k is a counter for fold number
            % 5s at 100Hz, is 500 data points. 80% for training = 400, 20%
            % validation = 100 data points
            RightTrain = zeros(400,8);
            LeftTrain = zeros(400,8);
            RightVal = zeros(100,8);
            LeftVal = zeros(100,8);

            for temp=1:numEpochs
               %Choose the start of the 100 data points randomly
               StartVal = randi([2,499]);
               if StartVal<401
                   RightEpochVal = ProcessedRightSide((temp-1)*500+StartVal:(temp-1)*500+StartVal+99,:);
                   LeftEpochVal = ProcessedLeftSide((temp-1)*500+StartVal:(temp-1)*500+StartVal+99,:);
                   RightEpochTr = [ProcessedRightSide((temp-1)*500+1:(temp-1)*500+StartVal-1,:); ProcessedRightSide((temp-1)*500+StartVal+100:temp*500,:)];
                   LeftEpochTr = [ProcessedLeftSide((temp-1)*500+1:(temp-1)*500+StartVal-1,:); ProcessedLeftSide((temp-1)*500+StartVal+100:temp*500,:)];
               else
                   RightEpochTr = ProcessedRightSide((temp-1)*500+StartVal-400:(temp-1)*500+StartVal-1,:);
                   LeftEpochTr = ProcessedLeftSide((temp-1)*500+StartVal-400:(temp-1)*500+StartVal-1,:);
                   RightEpochVal = [ProcessedRightSide((temp-1)*500+1:(temp-1)*500+StartVal-401,:); ProcessedRightSide((temp-1)*500+StartVal:temp*500,:)];
                   LeftEpochVal = [ProcessedLeftSide((temp-1)*500+1:(temp-1)*500+StartVal-401,:); ProcessedLeftSide((temp-1)*500+StartVal:temp*500,:)];
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

            save(strcat('Weight_Subj_', SubjID, '_Right_Tr_Fold_', num2str(k), '.mat' ), 'RightTrain')
            save(strcat('Weight_Subj_', SubjID, '_Left_Tr_Fold_', num2str(k), '.mat' ), 'LeftTrain')
            save(strcat('Weight_Subj_', SubjID, '_Right_Val_Fold_', num2str(k), '.mat' ), 'RightVal')
            save(strcat('Weight_Subj_', SubjID, '_Left_Val_Fold_', num2str(k), '.mat' ), 'LeftVal')

        end
    
    end
end