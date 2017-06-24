% to train Non-Negative Matrix Factorization and crossvalidate using the 20 folds for each subject
% Outcome will be a .mat file (NNMF_Subj_01.mat)

% 20160323 Written by Navid Shirzad

function NNMF_BatchProcess(SubjectIDs, NumberFolds, HandDominance)
    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end,

        DOF = 8;
        RightVAF = zeros(DOF,NumberFolds); LeftVAF = zeros(DOF,NumberFolds); %VAF for each number of synergies
        DeltaRightVAF = zeros(DOF,NumberFolds); DeltaLeftVAF = zeros(DOF,NumberFolds); %how much VAF changes when a synergy vector is added
        RightDOF_VAF = zeros(DOF,DOF,NumberFolds); LeftDOF_VAF = zeros(DOF,DOF,NumberFolds); %how VAF of each DOF changes as more synergies are added for each fold. first DOF is the synergies included, the second is the DOF being observed, third is the fold number
        ndimR = zeros(NumberFolds,1); ndimL = zeros(NumberFolds,1); %number of synergies (based on VAF requirements) for Tr data
        CAFRightTr = zeros(NumberFolds,1); CAFLeftTr = zeros(NumberFolds,1);
        AvgReconsErrRightVal = zeros(NumberFolds,1); AvgReconsErrLeftVal = zeros(NumberFolds,1); 
        AvgReconsErrRightTr = zeros(NumberFolds,1); AvgReconsErrLeftTr = zeros(NumberFolds,1);
        RegCoeffRightVal = zeros(NumberFolds,1);%Average correlation coefficient
        RegCoeffLeftVal = zeros(NumberFolds,1);
        RegCoeffRightTr = zeros(NumberFolds,1);
        RegCoeffLeftTr = zeros(NumberFolds,1);
        nCommonR = zeros(NumberFolds,1); nCommonL = zeros(NumberFolds,1);
        AllSynRTr = zeros(DOF,DOF,NumberFolds);
        AllSynLTr = zeros(DOF,DOF,NumberFolds);
        AllSynRVal = zeros(DOF,DOF,NumberFolds);
        AllSynLVal = zeros(DOF,DOF,NumberFolds);
        PCsRightTrAll = zeros(DOF,DOF,NumberFolds);
        PCsLeftTrAll = zeros(DOF,DOF,NumberFolds);
        PCsRightValAll = zeros(DOF,DOF,NumberFolds);
        PCsLeftValAll = zeros(DOF,DOF,NumberFolds);
        
        tic
        %there are k=20 folds
        for fold = 1:NumberFolds
            fold
            CurrentDirectory = cd;
            CurrentDirectoryUp = strrep(CurrentDirectory,'NNMF',''); %we are in PCA, data is saved up one level.
            %load the training and validation mat files
            load(strcat(CurrentDirectoryUp,'Weight_Subj_', SubjID, '_Right_Tr_Fold_', num2str(fold), '.mat' )); % loads it as 'RightTrain'
            load(strcat(CurrentDirectoryUp,'Weight_Subj_', SubjID, '_Left_Tr_Fold_', num2str(fold), '.mat' )); % 'LeftTrain'
            load(strcat(CurrentDirectoryUp,'Weight_Subj_', SubjID, '_Right_Val_Fold_', num2str(fold), '.mat' )); % 'RightVal'
            load(strcat(CurrentDirectoryUp,'Weight_Subj_', SubjID, '_Left_Val_Fold_', num2str(fold), '.mat' )); % 'LeftVal'
            
            %make sure data is non-negative: EMG data has been preprocessed
            %to make sure it is non-negative (in fact normalized to each
            %channel's max and saved as a percentage). So, yes, check!
            
            %to play around to see what is important
            RightTrain = RightTrain(:,1:DOF);
            LeftTrain = LeftTrain(:,1:DOF);
            RightVal = RightVal(:,1:DOF);
            LeftVal = LeftVal(:,1:DOF);

            
            %%%%check training is properly done (and record all the data so you can
            %%%%generate plots)
            
            GoodTrainR = 0; GoodTrainL=0; GoodTrainRVal = 0; GoodTrainLVal = 0;
            numSynergy = 1;
            while numSynergy < DOF 
                %perform NNMF on data
                true = 0; count = 1;
                while true == 0 && count < DOF && GoodTrainR == 0
                    [ScoresRightTrTemp, SynergiesRightTrAll] = nnmf(RightTrain, numSynergy); 
                    if rank(SynergiesRightTrAll) == numSynergy
                        true = 1; %not underfitting or stuck in local minima
                    end
                    count = count + 1;
                    if count == DOF
                        SynergiesRightTrAll = zeros(size(SynergiesRightTrAll));
                    end
                end
                true = 0; count = 1;
                while true == 0 && count < DOF && GoodTrainL == 0
                    [ScoresLeftTrTemp, SynergiesLeftTrAll] = nnmf(LeftTrain, numSynergy); 
                    if rank(SynergiesLeftTrAll) == numSynergy
                        true = 1; %not underfitting or stuck in local minima
                    end
                    count = count + 1;
                    if count == DOF
                        SynergiesLeftTrAll = zeros(size(SynergiesLeftTrAll));
                    end
                end
                true = 0; count = 1;
                while true == 0 && count < DOF && GoodTrainRVal == 0
                    [ScoresRightValTemp, SynergiesRightValAll] = nnmf(RightVal, numSynergy); 
                    if rank(SynergiesRightValAll) == numSynergy
                        true = 1; %not underfitting or stuck in local minima
                    end
                    count = count + 1;
                    if count == DOF
                        SynergiesRightValAll = zeros(size(SynergiesRightValAll));
                    end
                end
                true = 0; count = 1;
                while true == 0 && count < DOF && GoodTrainLVal == 0
                    [ScoresLeftValTemp, SynergiesLeftValAll] = nnmf(LeftVal, numSynergy);
                    if rank(SynergiesLeftValAll) == numSynergy
                        true = 1; %not underfitting or stuck in local minima
                    end
                    count = count + 1;
                    if count == DOF
                        SynergiesLeftValAll = zeros(size(SynergiesLeftValAll));
                    end
                end

                if GoodTrainR == 0 || GoodTrainL == 0
                    RightTrApprox = ScoresRightTrTemp * SynergiesRightTrAll;
                    RightVAF(numSynergy,fold) = 100*(1 - (sum(sum((RightTrain - RightTrApprox).^2,2),1)) / (sum(sum((RightTrain).^2,2),1))); %1-SSE/SST
                    LeftTrApprox = ScoresLeftTrTemp * SynergiesLeftTrAll;
                    LeftVAF(numSynergy,fold) = 100*(1 - (sum(sum((LeftTrain - LeftTrApprox).^2,2),1)) / (sum(sum((LeftTrain).^2,2),1))); %1-SSE/SST
                    if numSynergy ~= 1
                        DeltaRightVAF(numSynergy,fold)=RightVAF(numSynergy,fold)-RightVAF(numSynergy-1,fold);
                        DeltaLeftVAF(numSynergy,fold)=LeftVAF(numSynergy,fold)-LeftVAF(numSynergy-1,fold);
                    else
                        DeltaRightVAF(numSynergy,fold)=RightVAF(numSynergy,fold);
                        DeltaLeftVAF(numSynergy,fold)=LeftVAF(numSynergy,fold);
                    end
                    RightDOF_VAF(numSynergy,:,fold) = 100*(1 - sum((RightTrain - RightTrApprox).^2,1) ./ sum((RightTrain).^2,1));
                    LeftDOF_VAF(numSynergy,:,fold) = 100*(1 - sum((LeftTrain - LeftTrApprox).^2,1) ./ sum((LeftTrain).^2,1));
                end

                if GoodTrainR==0 & RightVAF(numSynergy,fold)>90 & DeltaRightVAF(numSynergy,fold)<5 & RightDOF_VAF(numSynergy,:,fold)>50 
                    ndimR(fold,1) = numSynergy;
                    GoodTrainR = 1
                    PCsRightTr = SynergiesRightTrAll;
                    PCsRightTrAll(1:size(PCsRightTr,1),:,fold) = PCsRightTr;
                    ScoresRightTr = ScoresRightTrTemp;
                    AllSynRTr(1:ndimR(fold,1),:, fold) = SynergiesRightTrAll;
                end
                if GoodTrainL==0 & LeftVAF(numSynergy,fold)>90 & DeltaLeftVAF(numSynergy,fold)<5 & LeftDOF_VAF(numSynergy,:,fold)>50 
                    ndimL(fold,1) = numSynergy;
                    GoodTrainL = 1
                    PCsLeftTr = SynergiesLeftTrAll;
                    PCsLeftTrAll(1:size(PCsLeftTr,1),:,fold) = PCsLeftTr;
                    ScoresLeftTr = ScoresLeftTrTemp;
                    AllSynLTr(1:ndimL(fold,1),:, fold) = SynergiesLeftTrAll;
                end
                
                if GoodTrainRVal == 0 || GoodTrainLVal == 0
                    RightValApprox = ScoresRightValTemp * SynergiesRightValAll;
                    ValRightVAF(numSynergy,fold) = 100*(1 - (sum(sum((RightVal - RightValApprox).^2,2),1)) / (sum(sum((RightVal).^2,2),1))); %1-SSE/SST
                    LeftValApprox = ScoresLeftValTemp * SynergiesLeftValAll;
                    ValLeftVAF(numSynergy,fold) = 100*(1 - (sum(sum((LeftVal - LeftValApprox).^2,2),1)) / (sum(sum((LeftVal).^2,2),1))); %1-SSE/SST
                    if numSynergy ~= 1
                        ValDeltaRightVAF(numSynergy,fold)=ValRightVAF(numSynergy,fold)-ValRightVAF(numSynergy-1,fold);
                        ValDeltaLeftVAF(numSynergy,fold)=ValLeftVAF(numSynergy,fold)-ValLeftVAF(numSynergy-1,fold);
                    else
                        ValDeltaRightVAF(numSynergy,fold)=ValRightVAF(numSynergy,fold);
                        ValDeltaLeftVAF(numSynergy,fold)=ValLeftVAF(numSynergy,fold);
                    end
                    ValRightDOF_VAF(numSynergy,:,fold) = 100*(1 - sum((RightVal - RightValApprox).^2,1) ./ sum((RightVal).^2,1));
                    ValLeftDOF_VAF(numSynergy,:,fold) = 100*(1 - sum((LeftVal - LeftValApprox).^2,1) ./ sum((LeftVal).^2,1));
                end
                
                if GoodTrainRVal==0 & ValRightVAF(numSynergy,fold)>90 & ValDeltaRightVAF(numSynergy,fold)<5 & ValRightDOF_VAF(numSynergy,:,fold)>50 
                    ndimRVal(fold,1) = numSynergy;
                    GoodTrainRVal = 1
                    PCsRightVal = SynergiesRightValAll;
                    PCsRightValAll(1:size(PCsRightVal,1),:,fold) = PCsRightVal;
                    AllSynRVal(1:ndimRVal(fold,1),:, fold) = SynergiesRightValAll;
                end
                if GoodTrainLVal==0 & ValLeftVAF(numSynergy,fold)>90 & ValDeltaLeftVAF(numSynergy,fold)<5 & ValLeftDOF_VAF(numSynergy,:,fold)>50 
                    ndimLVal(fold,1) = numSynergy;
                    GoodTrainLVal = 1
                    PCsLeftVal = SynergiesLeftValAll;
                    PCsLeftValAll(1:size(PCsLeftVal,1),:,fold) = PCsLeftVal;
                    AllSynLVal(1:ndimLVal(fold,1),:, fold) = SynergiesLeftValAll;
                end
                
                numSynergy = numSynergy+1;
                if numSynergy == DOF
                    if GoodTrainR==0 || GoodTrainL==0 || GoodTrainRVal==0 || GoodTrainLVal==0
                        numSynergy = 1;
                    end
                end
                if GoodTrainR==1 && GoodTrainL==1 && GoodTrainRVal==1 && GoodTrainLVal==1
                        numSynergy = DOF; %terminate training
                end
            end          
           
            %%%%validation

            %Covariance Accounted For CAF
            CAFRightTrAll = cumsum(diag(cov(ScoresRightTr)))/trace(cov(ScoresRightTr));
            CAFLeftTrAll = cumsum(diag(cov(ScoresLeftTr)))/trace(cov(ScoresLeftTr));
            CAFRightTr(fold,1) = CAFRightTrAll(ndimR(fold,1),1);
            CAFLeftTr(fold,1) = CAFLeftTrAll(ndimL(fold,1),1);

            %calculate the reconstruction of the data 
            RightEvalApprox = (RightVal / PCsRightTr) * PCsRightTr;
            LeftEvalApprox = (LeftVal / PCsLeftTr) * PCsLeftTr;
            RightTrApprox = ScoresRightTr * PCsRightTr;
            LeftTrApprox = ScoresLeftTr * PCsLeftTr;

            %avg reconstruction error of each data point (degrees)
            AvgReconsErrRightVal(fold,1) = (1/size(RightVal,1))*sum((1/DOF*sum((RightVal-RightEvalApprox).^2,2)).^0.5); %10 DOF
            AvgReconsErrLeftVal(fold,1) = (1/size(LeftVal,1))*sum((1/DOF*sum((LeftVal-LeftEvalApprox).^2,2)).^0.5); %10 DOF
            AvgReconsErrRightTr(fold,1) = (1/size(RightTrain,1))*sum((1/DOF*sum((RightTrain-RightTrApprox).^2,2)).^0.5); %10 DOF
            AvgReconsErrLeftTr(fold,1) = (1/size(LeftTrain,1))*sum((1/DOF*sum((LeftTrain -LeftTrApprox).^2,2)).^0.5); %10 DOF

            %correlation coefficient (slope)

            for temp = 1:DOF %10 DOFs
               [AllregcoeffREval(fold,temp), bint, r, rint, stats] = regress(RightEvalApprox(:,temp), RightVal(:,temp)); %slope
        % %        AllRSQRightVal(fold,temp) = stats(1,1); %RSquered
               [AllregcoeffLEval(fold,temp), bint, r, rint, stats] = regress(LeftEvalApprox(:,temp), LeftVal(:,temp));
        % %        AllRSQLeftVal(fold,temp) = stats(1,1);
               [AllregcoeffRTrain(fold,temp), bint, r, rint, stats] = regress(RightTrApprox(:,temp), RightTrain(:,temp));
        % %        AllRSQRightTr(fold,temp) = stats(1,1);
               [AllregcoeffLTrain(fold,temp), bint, r, rint, stats] = regress(LeftTrApprox(:,temp), LeftTrain(:,temp));
        % %        AllRSQLeftTr(fold,temp) = stats(1,1);
            end
            RegCoeffRightVal(fold,1) = sum(AllregcoeffREval(fold,:),2)/DOF;%Average correlation coefficient
            RegCoeffLeftVal(fold,1) = sum(AllregcoeffLEval(fold,:),2)/DOF;
            RegCoeffRightTr(fold,1) = sum(AllregcoeffRTrain(fold,:),2)/DOF;
            RegCoeffLeftTr(fold,1) = sum(AllregcoeffLTrain(fold,:),2)/DOF;

            %normalize synergy vectors
%             clear PCsRightValNorm
%             clear PCsLeftValNorm
%             clear PCsRightTrNorm
%             clear PCsLeftTrNorm
%             for temp=1:ndimRVal(fold,1)
%                 PCsRightValNorm(temp,:) = PCsRightVal(temp,:)/norm(PCsRightVal(temp,:)); %normalize vectors
%             end
%             for temp=1:ndimLVal(fold,1)
%                 PCsLeftValNorm(temp,:) = PCsLeftVal(temp,:)/norm(PCsLeftVal(temp,:)); %normalize vectors
%             end
%             for temp=1:ndimR(fold,1)
%                 PCsRightTrNorm(temp,:) = PCsRightTr(temp,:)/norm(PCsRightTr(temp,:)); %normalize vectors
%             end
%             for temp=1:ndimL(fold,1)
%                 PCsLeftTrNorm(temp,:) = PCsLeftTr(temp,:)/norm(PCsLeftTr(temp,:)); %normalize vectors
%             end

            %dim of common subspace
            for i=1:ndimLVal(fold,1)
                for j=1:ndimL(fold,1)
                    %AllLeftCommon(i,j) = 1 - (sum(sum((PCsLeftValNorm(:,i) - PCsLeftTrNorm(:,j)).^2,2),1)) / (sum(sum((PCsLeftTrNorm(:,i)).^2,2),1)); %1-SSE/SST
%                     [coeff, bint, r, rint, stats] = regress(PCsLeftValNorm(i,:)', PCsLeftTrNorm(j,:)');
%                     [coeff, bint, r, rint, stats] = regress(PCsLeftVal(i,:)', PCsLeftTr(j,:)');
%                     AllLeftCommon(i,j) = stats(1,1);
                    AllLeftCommonCoeff(i,j) = regress(PCsLeftVal(i,:)', PCsLeftTr(j,:)');
                end
            end
            for i=1:ndimRVal(fold,1)
                for j=1:ndimR(fold,1)
                    %AllRightCommon(i,j) = 1 - (sum(sum((PCsRightValNorm(:,i) - PCsRightTrNorm(:,j)).^2,2),1)) / (sum(sum((PCsRightValNorm(:,i)).^2,2),1)); %1-SSE/SST
%                     [coeff, bint, r, rint, stats] = regress(PCsRightValNorm(i,:)', PCsRightTrNorm(j,:)');
%                     [coeff, bint, r, rint, stats] = regress(PCsRightVal(i,:)', PCsRightTr(j,:)');
%                     AllRightCommon(i,j) = stats(1,1);
                    AllRightCommonCoeff(i,j) = regress(PCsRightVal(i,:)', PCsRightTr(j,:)');
                end
            end

            nCommonR(fold,1)=0; nCommonL(fold,1)=0;
            for i=1:ndimLVal(fold,1)
                if max(abs(AllLeftCommonCoeff(i,:)))>0.8
                    nCommonL(fold,1)=nCommonL(fold,1)+1;
                end
            end
            nCommonL(fold,1)=nCommonL(fold,1)/ndimLVal(fold,1);
            for i=1:ndimRVal(fold,1)
                if max(abs(AllRightCommonCoeff(i,:)))>0.8
                    nCommonR(fold,1)=nCommonR(fold,1)+1;
                end
            end
            nCommonR(fold,1)=nCommonR(fold,1)/ndimRVal(fold,1);
            
            for i=1:ndimL(fold,1)
                for j=1:ndimR(fold,1)
                    %AllLeftCommon(i,j) = 1 - (sum(sum((PCsLeftValNorm(:,i) - PCsLeftTrNorm(:,j)).^2,2),1)) / (sum(sum((PCsLeftTrNorm(:,i)).^2,2),1)); %1-SSE/SST
%                     [coeff, bint, r, rint, stats] = regress(PCsLeftTrNorm(:,i), PCsRightTrNorm(:,j));
                    BothCommonCoeff(i,j) = regress(PCsLeftTr(i,:)', PCsRightTr(j,:)');
                end
            end
            
            nCommon(fold,1)=0;
            for i=1:ndimL(fold,1)
                if max(BothCommonCoeff(i,:))>0.8
                    nCommon(fold,1)=nCommon(fold,1)+1;
                end
            end
            nCommon(fold,1)=nCommon(fold,1)/ndimLVal(fold,1);
        end
        timeElapsed = toc; 

        %%%find the average PC vectors
        for fold = 1:NumberFolds
            %do it for the first 5 PC vectors
            LeftPC1(fold,:) = AllSynLTr(1,:,fold);
            LeftPC2(fold,:) = AllSynLTr(2,:,fold);
            LeftPC3(fold,:) = AllSynLTr(3,:,fold);
            LeftPC4(fold,:) = AllSynLTr(4,:,fold);
            LeftPC5(fold,:) = AllSynLTr(5,:,fold);
            RightPC1(fold,:) = AllSynRTr(1,:,fold);
            RightPC2(fold,:) = AllSynRTr(2,:,fold);
            RightPC3(fold,:) = AllSynRTr(3,:,fold);
            RightPC4(fold,:) = AllSynRTr(4,:,fold);
            RightPC5(fold,:) = AllSynRTr(5,:,fold);           
        end


        LeftPC1Avg = zeros(1, DOF);
        LeftPC2Avg = zeros(1, DOF);
        LeftPC3Avg = zeros(1, DOF);
        LeftPC4Avg = zeros(1, DOF);
        LeftPC5Avg = zeros(1, DOF);
        RightPC1Avg = zeros(1, DOF);
        RightPC2Avg = zeros(1, DOF);
        RightPC3Avg = zeros(1, DOF);
        RightPC4Avg = zeros(1, DOF);
        RightPC5Avg = zeros(1, DOF);

        for fold = 1:NumberFolds
            %do it for the first 5 PC vectors
            LeftPC1Avg= LeftPC1Avg+LeftPC1(fold,:);
            LeftPC2Avg= LeftPC2Avg+LeftPC2(fold,:);
            LeftPC3Avg= LeftPC3Avg+LeftPC3(fold,:);
            LeftPC4Avg= LeftPC4Avg+LeftPC4(fold,:);
            LeftPC5Avg= LeftPC5Avg+LeftPC5(fold,:);
            RightPC1Avg= RightPC1Avg+RightPC1(fold,:);
            RightPC2Avg= RightPC2Avg+RightPC2(fold,:);
            RightPC3Avg= RightPC3Avg+RightPC3(fold,:);
            RightPC4Avg= RightPC4Avg+RightPC4(fold,:);
            RightPC5Avg= RightPC5Avg+RightPC5(fold,:);

        end

        LeftPC1Avg= LeftPC1Avg/NumberFolds;
        LeftPC2Avg= LeftPC2Avg/NumberFolds;
        LeftPC3Avg= LeftPC3Avg/NumberFolds;
        LeftPC4Avg= LeftPC4Avg/NumberFolds;
        LeftPC5Avg= LeftPC5Avg/NumberFolds;
        RightPC1Avg= RightPC1Avg/NumberFolds;
        RightPC2Avg= RightPC2Avg/NumberFolds;
        RightPC3Avg= RightPC3Avg/NumberFolds;
        RightPC4Avg= RightPC4Avg/NumberFolds;
        RightPC5Avg= RightPC5Avg/NumberFolds;

               
        %% save the analysis outcome
            if HandDominance(subjectcounter) == 1
                SubjDomHand = 'R';
            else
                SubjDomHand = 'L';
            end
            
                     
            save(strcat('NNMF_', 'Weight_Subj_', SubjID, '.mat' ), 'RightVAF', 'LeftVAF', ...
            'DeltaRightVAF', 'DeltaLeftVAF', 'RightDOF_VAF', 'DeltaLeftVAF', ...
            'ndimR', 'ndimL', 'CAFRightTr', 'CAFLeftTr', ...
            'AvgReconsErrRightVal', 'AvgReconsErrLeftVal', 'AvgReconsErrRightTr', 'AvgReconsErrLeftTr', ...
            'RegCoeffRightVal', 'RegCoeffLeftVal', 'RegCoeffRightTr', 'RegCoeffLeftTr', ...
            'nCommonR', 'nCommonL', 'nCommon',  ...
            'RightPC1Avg', 'RightPC2Avg', 'RightPC3Avg', 'RightPC4Avg', 'RightPC5Avg', ...
            'LeftPC1Avg', 'LeftPC2Avg', 'LeftPC3Avg', 'LeftPC4Avg', 'LeftPC5Avg', ...
            'SubjDomHand', 'timeElapsed', ...
            'PCsRightTrAll', 'PCsLeftTrAll', 'PCsRightValAll', 'PCsLeftValAll');
        
        
        %% Plotting function(when you call the main function, edit to define which hand is dominant/stronger)
        
        %         figure()
%         subplot(1,2,1)
%         bar(LeftPC1,'DisplayName','LeftPC1')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Left Arm')
%         ylabel('Normalized Weight of Each DOF')
%         subplot(1,2,2)
%         bar(RightPC1,'DisplayName','RightPC1')
%         xlabel('Muscle DOFs of Right Arm')
%         axis([0 9 0 1])
%         title(strcat('Subj_', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 1                                                           .'))
%         
%         figure()
%         subplot(1,2,1)
%         bar(LeftPC2,'DisplayName','LeftPC2')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Left Arm')
%         ylabel('Normalized Weight of Each DOF')
%         subplot(1,2,2)
%         bar(RightPC2,'DisplayName','RightPC2')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Right Arm')
%         title(strcat('Subj_', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 2                                                           .'))
%         
%         figure()
%         subplot(1,2,1)
%         bar(LeftPC3,'DisplayName','LeftPC3')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Left Arm')
%         ylabel('Normalized Weight of Each DOF')
%         subplot(1,2,2)
%         bar(RightPC3,'DisplayName','RightPC3')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Right Arm')
%         title(strcat('Subj_', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 3                                                           .'))
%         
%         figure()
%         subplot(1,2,1)
%         bar(LeftPC4,'DisplayName','LeftPC4')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Left Arm')
%         ylabel('Normalized Weight of Each DOF')
%         subplot(1,2,2)
%         bar(RightPC4,'DisplayName','RightPC4')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Right Arm')
%         title(strcat('Subj_', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 4                                                           .'))
%         
%         figure()
%         subplot(1,2,1)
%         bar(LeftPC5,'DisplayName','LeftPC5')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Left Arm')
%         ylabel('Normalized Weight of Each DOF')
%         subplot(1,2,2)
%         bar(RightPC5,'DisplayName','RightPC5')
%         axis([0 9 0 1])
%         xlabel('Muscle DOFs of Right Arm')
%         title(strcat('Subj_', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 5                                                           .'))
        
        figure()
        subplot(5,2,1)
        bar(LeftPC1','DisplayName','LeftPC1')
        axis([0 DOF+1 -1 1])
        subplot(5,2,2)
        bar(RightPC1','DisplayName','RightPC1')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 1                                                           .'))
        
        subplot(5,2,3)
        bar(LeftPC2','DisplayName','LeftPC2')
        axis([0 DOF+1 -1 1])
        subplot(5,2,4)
        bar(RightPC2','DisplayName','RightPC2')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 2                                                           .'))
        
        
        subplot(5,2,5)
        bar(LeftPC3','DisplayName','LeftPC3')
        axis([0 DOF+1 -1 1])
        ylabel('Normalized Weight of Each DOF')
        subplot(5,2,6)
        bar(RightPC3','DisplayName','RightPC3')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 3                                                           .'))
        
        subplot(5,2,7)
        bar(LeftPC4','DisplayName','LeftPC4')
        axis([0 DOF+1 -1 1])
        subplot(5,2,8)
        bar(RightPC4','DisplayName','RightPC4')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 4                                                           .'))
        
        subplot(5,2,9)
        bar(LeftPC5','DisplayName','LeftPC5')
        axis([0 DOF+1 -1 1])
        xlabel('Muscle DOFs of Left Arm')
        subplot(5,2,10)
        bar(RightPC5','DisplayName','RightPC5')
        axis([0 DOF+1 -1 1])
        xlabel('Muscle DOFs of Right Arm')
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 5                                                           .'))

         figure()
        subplot(5,2,1)
        bar(LeftPC1Avg,'DisplayName','LeftPC1')
        axis([0 DOF+1 -1 1])
        subplot(5,2,2)
        bar(RightPC1Avg,'DisplayName','RightPC1')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 1                                                           .'))
        
        subplot(5,2,3)
        bar(LeftPC2Avg,'DisplayName','LeftPC2')
        axis([0 DOF+1 -1 1])
        subplot(5,2,4)
        bar(RightPC2Avg,'DisplayName','RightPC2')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 2                                                           .'))
        
        
        subplot(5,2,5)
        bar(LeftPC3Avg,'DisplayName','LeftPC3')
        axis([0 DOF+1 -1 1])
        ylabel('Normalized Weight of Each DOF')
        subplot(5,2,6)
        bar(RightPC3Avg,'DisplayName','RightPC3')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 3                                                           .'))
        
        subplot(5,2,7)
        bar(LeftPC4Avg,'DisplayName','LeftPC4')
        axis([0 DOF+1 -1 1])
        subplot(5,2,8)
        bar(RightPC4Avg,'DisplayName','RightPC4')
        axis([0 DOF+1 -1 1])
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 4                                                           .'))
        
        subplot(5,2,9)
        bar(LeftPC5Avg,'DisplayName','LeftPC5')
        axis([0 DOF+1 -1 1])
        xlabel('Muscle DOFs of Left Arm')
        subplot(5,2,10)
        bar(RightPC5Avg,'DisplayName','RightPC5')
        axis([0 DOF+1 -1 1])
        xlabel('Muscle DOFs of Right Arm')
        title(strcat('Subj ', num2str(SubjectIDs(subjectcounter)),' Muscle Synergy Vector 5                                                           .'))
    
    end
end