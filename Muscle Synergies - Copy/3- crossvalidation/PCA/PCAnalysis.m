% to train PCA and crossvalidate using the 10 folds for each subject
%change subjectID

% 20150911 Written by Navid Shirzad

clear; clc;

foldNumber = 20;
subjectID = 10; 
DOF = 10;
RightVAF = zeros(DOF,foldNumber); LeftVAF = zeros(DOF,foldNumber); %VAF for each number of synergies
DeltaRightVAF = zeros(DOF,foldNumber); DeltaLeftVAF = zeros(DOF,foldNumber); %how much VAF changes when a synergy vector is added
RightDOF_VAF = zeros(DOF,DOF,foldNumber); LeftDOF_VAF = zeros(DOF,DOF,foldNumber); %how VAF of each DOF changes as more synergies are added for each fold. first DOF is the synergies included, the second is the DOF being observed, third is the fold number
ndimR = zeros(foldNumber,1); ndimL = zeros(foldNumber,1);
CAFRightTr = zeros(foldNumber,1); CAFLeftTr = zeros(foldNumber,1);
AvgReconsErrRightVal = zeros(foldNumber,1); AvgReconsErrLeftVal = zeros(foldNumber,1); AvgReconsErrRightTr = zeros(foldNumber,1); AvgReconsErrLeftTr = zeros(foldNumber,1);

%there are k=20 folds
for fold = 1:foldNumber
	%load the training and validation mat files
	load(strcat('Y', num2str(subjectID), '_Right_Tr_Fold_', num2str(fold), '.mat' )); % loads it as 'RightTrain'
    load(strcat('Y', num2str(subjectID), '_Left_Tr_Fold_', num2str(fold), '.mat' )); % 'LeftTrain'
    load(strcat('Y', num2str(subjectID), '_Right_Val_Fold_', num2str(fold), '.mat' )); % 'RightVal'
    load(strcat('Y', num2str(subjectID), '_Left_Val_Fold_', num2str(fold), '.mat' )); % 'LeftVal'
      
    %perform PCA on training data
    [PCsRightTrAll(:,:,fold),ScoresRightTr,VarRightTr] = pca(RightTrain,'Economy',true,'Centered',false); %false for economy returns all pcs. no dimensionality reduction
    [PCsLeftTrAll(:,:,fold),ScoresLeftTr,VarLeftTr] = pca(LeftTrain,'Economy',true,'Centered',false); %false for economy returns all pcs. no dimensionality reduction. Centered false does not subtract column means.
    [PCsRightValAll(:,:,fold),ScoresRightVal,VarRightVal] = pca(RightVal,'Economy',true,'Centered',false); %false for economy returns all pcs. no dimensionality reduction
    [PCsLeftValAll(:,:,fold),ScoresLeftVal,VarLeftVal] = pca(LeftVal,'Economy',true,'Centered',false);
    
    PCsRightTr = PCsRightTrAll(:,:,fold);
    PCsLeftTr = PCsLeftTrAll(:,:,fold);
    PCsRightVal = PCsRightValAll(:,:,fold);
    PCsLeftVal = PCsLeftValAll(:,:,fold);
    %%%%check training is properly done (and record all the data so you can
    %%%%generate plots)
    GoodTrainR = 0; GoodTrainL=0;
    for numSynergy = 1:DOF 
        RightTrApprox = ScoresRightTr(:,1:numSynergy)*PCsRightTr(:,1:numSynergy)';
        RightVAF(numSynergy,fold) = 100*(1 - (sum(sum((RightTrain - RightTrApprox).^2,2),1)) / (sum(sum((RightTrain).^2,2),1))); %1-SSE/SST
        LeftTrApprox = ScoresLeftTr(:,1:numSynergy)*PCsLeftTr(:,1:numSynergy)';
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
        
        if RightVAF(numSynergy,fold)>90 & DeltaRightVAF(numSynergy,fold)<3 & RightDOF_VAF(numSynergy,4:DOF,fold)>50 & GoodTrainR==0
            ndimR(fold,1) = numSynergy;
            GoodTrainR = 1;
        end
        if LeftVAF(numSynergy,fold)>90 & DeltaLeftVAF(numSynergy,fold)<3 & LeftDOF_VAF(numSynergy,4:DOF,fold)>50 & GoodTrainL==0
            ndimL(fold,1) = numSynergy;
            GoodTrainL = 1;
        end        
    end
    
    %%%%validation
    
    %Covariance Accounted For CAF
    CAFRightTrAll(1,:) = cumsum(diag(cov(ScoresRightTr)))/trace(cov(ScoresRightTr));
    CAFLeftTrAll(1,:) = cumsum(diag(cov(ScoresLeftTr)))/trace(cov(ScoresLeftTr));
    CAFRightTr(fold,1) = CAFRightTrAll(1,ndimR(fold,1));
    CAFLeftTr(fold,1) = CAFLeftTrAll(1,ndimL(fold,1));
    
    %demean Eval data and Training data
    meansR=mean(RightVal);
    meansL=mean(LeftVal);
    meansMatrixREval=[ones(size(RightVal,1),1)*meansR(1,1) ones(size(RightVal,1),1)*meansR(1,2) ones(size(RightVal,1),1)*meansR(1,3) ...
        ones(size(RightVal,1),1)*meansR(1,4) ones(size(RightVal,1),1)*meansR(1,5) ones(size(RightVal,1),1)*meansR(1,6) ...
        ones(size(RightVal,1),1)*meansR(1,7) ones(size(RightVal,1),1)*meansR(1,8) ones(size(RightVal,1),1)*meansR(1,9) ones(size(RightVal,1),1)*meansR(1,10)];
    meansMatrixLEval=[ones(size(LeftVal,1),1)*meansL(1,1) ones(size(LeftVal,1),1)*meansL(1,2) ones(size(LeftVal,1),1)*meansL(1,3) ...
        ones(size(LeftVal,1),1)*meansL(1,4) ones(size(LeftVal,1),1)*meansL(1,5) ones(size(LeftVal,1),1)*meansL(1,6) ...
        ones(size(LeftVal,1),1)*meansL(1,7) ones(size(LeftVal,1),1)*meansL(1,8) ones(size(LeftVal,1),1)*meansL(1,9) ones(size(LeftVal,1),1)*meansL(1,10)];

%     meansR=mean(RightTrain);
%     meansL=mean(LeftTrain);
%     meansMatrixRTrain=[ones(size(RightTrain,1),1)*meansR(1,1) ones(size(RightTrain,1),1)*meansR(1,2) ones(size(RightTrain,1),1)*meansR(1,3) ...
%         ones(size(RightTrain,1),1)*meansR(1,4) ones(size(RightTrain,1),1)*meansR(1,5) ones(size(RightTrain,1),1)*meansR(1,6) ...
%         ones(size(RightTrain,1),1)*meansR(1,7) ones(size(RightTrain,1),1)*meansR(1,8) ones(size(RightTrain,1),1)*meansR(1,9) ones(size(RightTrain,1),1)*meansR(1,10)];
%     meansMatrixLTrain=[ones(size(LeftTrain,1),1)*meansL(1,1) ones(size(LeftTrain,1),1)*meansL(1,2) ones(size(LeftTrain,1),1)*meansL(1,3) ...
%         ones(size(LeftTrain,1),1)*meansL(1,4) ones(size(LeftTrain,1),1)*meansL(1,5) ones(size(LeftTrain,1),1)*meansL(1,6) ...
%         ones(size(LeftTrain,1),1)*meansL(1,7) ones(size(LeftTrain,1),1)*meansL(1,8) ones(size(LeftTrain,1),1)*meansL(1,9) ones(size(LeftTrain,1),1)*meansL(1,10)];
    
    %calculate the reconstruction of the data 
    RightEvalApprox = (RightVal - meansMatrixREval) * [PCsRightTr(:,1:ndimR(fold,1)) zeros(size(PCsRightTr,1),10-ndimR(fold,1))] * PCsRightTr';
    LeftEvalApprox = (LeftVal - meansMatrixLEval) * [PCsLeftTr(:,1:ndimL(fold,1)) zeros(size(PCsLeftTr,1),10-ndimL(fold,1))] * PCsLeftTr';
    RightTrApprox = ScoresRightTr(:,1:ndimR(fold,1))*PCsRightTr(:,1:ndimR(fold,1))';
    LeftTrApprox = ScoresLeftTr(:,1:ndimL(fold,1))*PCsLeftTr(:,1:ndimL(fold,1))';
    
    %avg reconstruction error of each data point (degrees)
    AvgReconsErrRightVal(fold,1) = (1/size(RightVal,1))*sum((1/DOF*sum(((RightVal - meansMatrixREval)-RightEvalApprox).^2,2)).^0.5); %10 DOF
    AvgReconsErrLeftVal(fold,1) = (1/size(LeftVal,1))*sum((1/DOF*sum(((LeftVal - meansMatrixLEval)-LeftEvalApprox).^2,2)).^0.5); %10 DOF
    AvgReconsErrRightTr(fold,1) = (1/size(RightTrain,1))*sum((1/DOF*sum((RightTrain-RightTrApprox).^2,2)).^0.5); %10 DOF
    AvgReconsErrLeftTr(fold,1) = (1/size(LeftTrain,1))*sum((1/DOF*sum((LeftTrain -LeftTrApprox).^2,2)).^0.5); %10 DOF
    
    %correlation coefficient (slope)
    
    for temp = 1:DOF %10 DOFs
       [AllregcoeffREval(fold,temp), bint, r, rint, stats] = regress(RightEvalApprox(:,temp), RightVal(:,temp) - meansMatrixREval(:,temp)); %slope
% %        AllRSQRightVal(fold,temp) = stats(1,1); %RSquered
       [AllregcoeffLEval(fold,temp), bint, r, rint, stats] = regress(LeftEvalApprox(:,temp), LeftVal(:,temp) - meansMatrixLEval(:,temp));
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
    
    
    %find the number of synergies for the Val dataset
    GoodTrainR = 0; GoodTrainL=0;numSynergy=1;
    while  GoodTrainR == 0 | GoodTrainL==0;
        RightValApprox = ScoresRightVal(:,1:numSynergy)*PCsRightVal(:,1:numSynergy)';
        ValRightVAF(numSynergy,fold) = 100*(1 - (sum(sum((RightVal - RightValApprox).^2,2),1)) / (sum(sum((RightVal).^2,2),1))); %1-SSE/SST
        LeftValApprox = ScoresLeftVal(:,1:numSynergy)*PCsLeftVal(:,1:numSynergy)';
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
        
        if ValRightVAF(numSynergy,fold)>90 & ValDeltaRightVAF(numSynergy,fold)<3 & ValRightDOF_VAF(numSynergy,4:DOF,fold)>50 & GoodTrainR==0
            ndimRVal(fold,1) = numSynergy;
            GoodTrainR = 1;
        end
        if ValLeftVAF(numSynergy,fold)>90 & ValDeltaLeftVAF(numSynergy,fold)<3 & ValLeftDOF_VAF(numSynergy,4:DOF,fold)>50 & GoodTrainL==0
            ndimLVal(fold,1) = numSynergy;
            GoodTrainL = 1;
        end 
        numSynergy=numSynergy+1;
    end
    
    %normalize synergy vectors
    for temp=1:DOF
        PCsRightValNorm(:,temp) = PCsRightVal(:,temp)/norm(PCsRightVal(:,temp)); %normalize vectors
        PCsLeftValNorm(:,temp) = PCsLeftVal(:,temp)/norm(PCsLeftVal(:,temp)); %normalize vectors
        PCsRightTrNorm(:,temp) = PCsRightTr(:,temp)/norm(PCsRightTr(:,temp)); %normalize vectors
        PCsLeftTrNorm(:,temp) = PCsLeftTr(:,temp)/norm(PCsLeftTr(:,temp)); %normalize vectors
    end
    
    %dim of common subspace
    for i=1:ndimLVal(fold,1)
        for j=1:ndimL(fold,1)
            %AllLeftCommon(i,j) = 1 - (sum(sum((PCsLeftValNorm(:,i) - PCsLeftTrNorm(:,j)).^2,2),1)) / (sum(sum((PCsLeftTrNorm(:,i)).^2,2),1)); %1-SSE/SST
            [coeff, bint, r, rint, stats] = regress(PCsLeftValNorm(:,i), PCsLeftTrNorm(:,j));
            AllLeftCommon(i,j) = stats(1,1);
        end
    end
    for i=1:ndimRVal(fold,1)
        for j=1:ndimR(fold,1)
            %AllRightCommon(i,j) = 1 - (sum(sum((PCsRightValNorm(:,i) - PCsRightTrNorm(:,j)).^2,2),1)) / (sum(sum((PCsRightValNorm(:,i)).^2,2),1)); %1-SSE/SST
            [coeff, bint, r, rint, stats] = regress(PCsRightValNorm(:,i), PCsRightTrNorm(:,j));
            AllRightCommon(i,j) = stats(1,1);
        end
    end
    
    nCommonR(fold,1)=0; nCommonL(fold,1)=0;
    for i=1:ndimLVal(fold,1)
        if max(AllLeftCommon(i,:))>0.9
            nCommonL(fold,1)=nCommonL(fold,1)+1;
        end
    end
    nCommonL(fold,1)=nCommonL(fold,1)/ndimLVal(fold,1);
    for i=1:ndimRVal(fold,1)
        if max(AllRightCommon(i,:))>0.9
            nCommonR(fold,1)=nCommonR(fold,1)+1;
        end
    end
    nCommonR(fold,1)=nCommonR(fold,1)/ndimRVal(fold,1);

end

    
%%%find the average PC vectors
for fold = 1:foldNumber
    %do it for the first 5 PC vectors
    
    if size(find(PCsLeftTrAll(:,1,fold)<0),1)<5
        LeftPC1(:,fold)= PCsLeftTrAll(:,1,fold)/norm(PCsLeftTrAll(:,1,fold));
    else
        LeftPC1(:,fold)= -PCsLeftTrAll(:,1,fold)/norm(PCsLeftTrAll(:,1,fold));
    end
    if size(find(PCsLeftTrAll(:,2,fold)<0),1)<5
        LeftPC2(:,fold)= PCsLeftTrAll(:,2,fold)/norm(PCsLeftTrAll(:,2,fold));
    else
        LeftPC2(:,fold)= -PCsLeftTrAll(:,2,fold)/norm(PCsLeftTrAll(:,2,fold));
    end
    if size(find(PCsLeftTrAll(:,3,fold)<0),1)<5
        LeftPC3(:,fold)= PCsLeftTrAll(:,3,fold)/norm(PCsLeftTrAll(:,3,fold));
    else
        LeftPC3(:,fold)= -PCsLeftTrAll(:,3,fold)/norm(PCsLeftTrAll(:,3,fold));
    end
    if size(find(PCsLeftTrAll(:,4,fold)<0),1)<5
        LeftPC4(:,fold)= PCsLeftTrAll(:,4,fold)/norm(PCsLeftTrAll(:,4,fold));
    else
        LeftPC4(:,fold)= -PCsLeftTrAll(:,4,fold)/norm(PCsLeftTrAll(:,4,fold));
    end
    if size(find(PCsLeftTrAll(:,5,fold)<0),1)<5
        LeftPC5(:,fold)= PCsLeftTrAll(:,5,fold)/norm(PCsLeftTrAll(:,5,fold));
    else
        LeftPC5(:,fold)= -PCsLeftTrAll(:,5,fold)/norm(PCsLeftTrAll(:,5,fold));
    end
    
    if size(find(PCsRightTrAll(:,1,fold)<0),1)<5
        RightPC1(:,fold)= PCsRightTrAll(:,1,fold)/norm(PCsRightTrAll(:,1,fold));
    else
        RightPC1(:,fold)= -PCsRightTrAll(:,1,fold)/norm(PCsRightTrAll(:,1,fold));
    end
    if size(find(PCsRightTrAll(:,2,fold)<0),1)<5
        RightPC2(:,fold)= PCsRightTrAll(:,2,fold)/norm(PCsRightTrAll(:,2,fold));
    else
        RightPC2(:,fold)= -PCsRightTrAll(:,2,fold)/norm(PCsRightTrAll(:,2,fold));
    end
    if size(find(PCsRightTrAll(:,3,fold)<0),1)<5
        RightPC3(:,fold)= PCsRightTrAll(:,3,fold)/norm(PCsRightTrAll(:,3,fold));
    else
        RightPC3(:,fold)= -PCsRightTrAll(:,3,fold)/norm(PCsRightTrAll(:,3,fold));
    end
    if size(find(PCsRightTrAll(:,4,fold)<0),1)<5
        RightPC4(:,fold)= PCsRightTrAll(:,4,fold)/norm(PCsRightTrAll(:,4,fold));
    else
        RightPC4(:,fold)= -PCsRightTrAll(:,4,fold)/norm(PCsRightTrAll(:,4,fold));
    end
    if size(find(PCsRightTrAll(:,5,fold)<0),1)<5
        RightPC5(:,fold)= PCsRightTrAll(:,5,fold)/norm(PCsRightTrAll(:,5,fold));
    else
        RightPC5(:,fold)= -PCsRightTrAll(:,5,fold)/norm(PCsRightTrAll(:,5,fold));
    end
    
end


LeftPC1Avg(:,1)= zeros(DOF,1);
LeftPC2Avg(:,1)= zeros(DOF,1);
LeftPC3Avg(:,1)= zeros(DOF,1);
LeftPC4Avg(:,1)= zeros(DOF,1);
LeftPC5Avg(:,1)= zeros(DOF,1);
RightPC1Avg(:,1)= zeros(DOF,1);
RightPC2Avg(:,1)= zeros(DOF,1);
RightPC3Avg(:,1)= zeros(DOF,1);
RightPC4Avg(:,1)= zeros(DOF,1);
RightPC5Avg(:,1)= zeros(DOF,1);

for fold = 1:foldNumber
    %do it for the first 5 PC vectors
    LeftPC1Avg(:,1)= LeftPC1Avg(:,1)+LeftPC1(:,fold);
    LeftPC2Avg(:,1)= LeftPC2Avg(:,1)+LeftPC2(:,fold);
    LeftPC3Avg(:,1)= LeftPC3Avg(:,1)+LeftPC3(:,fold);
    LeftPC4Avg(:,1)= LeftPC4Avg(:,1)+LeftPC4(:,fold);
    LeftPC5Avg(:,1)= LeftPC5Avg(:,1)+LeftPC5(:,fold);
    RightPC1Avg(:,1)= RightPC1Avg(:,1)+RightPC1(:,fold);
    RightPC2Avg(:,1)= RightPC2Avg(:,1)+RightPC2(:,fold);
    RightPC3Avg(:,1)= RightPC3Avg(:,1)+RightPC3(:,fold);
    RightPC4Avg(:,1)= RightPC4Avg(:,1)+RightPC4(:,fold);
    RightPC5Avg(:,1)= RightPC5Avg(:,1)+RightPC5(:,fold);
    
end

LeftPC1Avg(:,1)= LeftPC1Avg(:,1)/foldNumber;
LeftPC2Avg(:,1)= LeftPC2Avg(:,1)/foldNumber;
LeftPC3Avg(:,1)= LeftPC3Avg(:,1)/foldNumber;
LeftPC4Avg(:,1)= LeftPC4Avg(:,1)/foldNumber;
LeftPC5Avg(:,1)= LeftPC5Avg(:,1)/foldNumber;
RightPC1Avg(:,1)= RightPC1Avg(:,1)/foldNumber;
RightPC2Avg(:,1)= RightPC2Avg(:,1)/foldNumber;
RightPC3Avg(:,1)= RightPC3Avg(:,1)/foldNumber;
RightPC4Avg(:,1)= RightPC4Avg(:,1)/foldNumber;
RightPC5Avg(:,1)= RightPC5Avg(:,1)/foldNumber;


if size(find(LeftPC1Avg(:,1)<0),1)>4
    LeftPC1Avg(:,1)= -LeftPC1Avg(:,1)/foldNumber;
end
if size(find(LeftPC2Avg(:,1)<0),1)>4
    LeftPC2Avg(:,1)= -LeftPC2Avg(:,1)/foldNumber;
end
if size(find(LeftPC3Avg(:,1)<0),1)>4
    LeftPC3Avg(:,1)= -LeftPC3Avg(:,1)/foldNumber;
end
if size(find(LeftPC4Avg(:,1)<0),1)>4
    LeftPC4Avg(:,1)= -LeftPC4Avg(:,1)/foldNumber;
end
if size(find(LeftPC5Avg(:,1)<0),1)>4
    LeftPC5Avg(:,1)= -LeftPC5Avg(:,1)/foldNumber;
end


for fold=1:foldNumber
   if dot(LeftPC1Avg(:,1),LeftPC1(:,fold))<0
       LeftPC1(:,fold) = -LeftPC1(:,fold);
   end
    if dot(LeftPC2Avg(:,1),LeftPC2(:,fold))<0
       LeftPC2(:,fold) = -LeftPC2(:,fold);
    end
   if dot(LeftPC3Avg(:,1),LeftPC3(:,fold))<0
       LeftPC3(:,fold) = -LeftPC3(:,fold);
   end
   if dot(LeftPC4Avg(:,1),LeftPC4(:,fold))<0
       LeftPC4(:,fold) = -LeftPC4(:,fold);
   end
   if dot(LeftPC5Avg(:,1),LeftPC5(:,fold))<0
       LeftPC5(:,fold) = -LeftPC5(:,fold);
   end
   
   if dot(LeftPC1Avg(:,1),RightPC1(:,fold))<0
       RightPC1(:,fold) = -RightPC1(:,fold);
   end
    if dot(LeftPC2Avg(:,1),RightPC2(:,fold))<0
       RightPC2(:,fold) = -RightPC2(:,fold);
    end
   if dot(LeftPC3Avg(:,1),RightPC3(:,fold))<0
       RightPC3(:,fold) = -RightPC3(:,fold);
   end
   if dot(LeftPC4Avg(:,1),RightPC4(:,fold))<0
       RightPC4(:,fold) = -RightPC4(:,fold);
   end
   if dot(LeftPC5Avg(:,1),RightPC5(:,fold))<0
       RightPC5(:,fold) = -RightPC5(:,fold);
   end
end



%%use when right hand is not affected

% % % % figure()
% % % % subplot(1,2,1)
% % % % bar(LeftPC1,'DisplayName','LeftPC1')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(1,2,2)
% % % % bar(RightPC1,'DisplayName','RightPC1')
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % axis([0 11 -1 1])
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 1                                                           .'))
% % % % 
% % % % figure()
% % % % subplot(1,2,1)
% % % % bar(LeftPC2,'DisplayName','LeftPC2')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(1,2,2)
% % % % bar(RightPC2,'DisplayName','RightPC2')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 2                                                           .'))
% % % % 
% % % % figure()
% % % % subplot(1,2,1)
% % % % bar(LeftPC3,'DisplayName','LeftPC3')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(1,2,2)
% % % % bar(RightPC3,'DisplayName','RightPC3')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 3                                                           .'))
% % % % 
% % % % figure()
% % % % subplot(1,2,1)
% % % % bar(LeftPC4,'DisplayName','LeftPC4')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(1,2,2)
% % % % bar(RightPC4,'DisplayName','RightPC4')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 4                                                           .'))
% % % % 
% % % % figure()
% % % % subplot(1,2,1)
% % % % bar(LeftPC5,'DisplayName','LeftPC5')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(1,2,2)
% % % % bar(RightPC5,'DisplayName','RightPC5')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 5                                                           .'))
% % % % 
% % % % figure()
% % % % subplot(5,2,1)
% % % % bar(LeftPC1,'DisplayName','LeftPC1')
% % % % axis([0 11 -1 1])
% % % % subplot(5,2,2)
% % % % bar(RightPC1,'DisplayName','RightPC1')
% % % % axis([0 11 -1 1])
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 1                                                           .'))
% % % % 
% % % % subplot(5,2,3)
% % % % bar(LeftPC2,'DisplayName','LeftPC2')
% % % % axis([0 11 -1 1])
% % % % subplot(5,2,4)
% % % % bar(RightPC2,'DisplayName','RightPC2')
% % % % axis([0 11 -1 1])
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 2                                                           .'))
% % % % 
% % % % 
% % % % subplot(5,2,5)
% % % % bar(LeftPC3,'DisplayName','LeftPC3')
% % % % axis([0 11 -1 1])
% % % % ylabel('Normalized Weight of Each DOF')
% % % % subplot(5,2,6)
% % % % bar(RightPC3,'DisplayName','RightPC3')
% % % % axis([0 11 -1 1])
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 3                                                           .'))
% % % % 
% % % % subplot(5,2,7)
% % % % bar(LeftPC4,'DisplayName','LeftPC4')
% % % % axis([0 11 -1 1])
% % % % subplot(5,2,8)
% % % % bar(RightPC4,'DisplayName','RightPC4')
% % % % axis([0 11 -1 1])
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 4                                                           .'))
% % % % 
% % % % subplot(5,2,9)
% % % % bar(LeftPC5,'DisplayName','LeftPC5')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Impaired Arm')
% % % % subplot(5,2,10)
% % % % bar(RightPC5,'DisplayName','RightPC5')
% % % % axis([0 11 -1 1])
% % % % xlabel('Kinematic DOFs of Unimpaired Arm')
% % % % title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 5                                                           .'))


%% use when left hand is not impaired
figure()
subplot(1,2,1)
bar(RightPC1,'DisplayName','RightPC1')
xlabel('Kinematic DOFs of Impaired Arm')
axis([0 11 -1 1])
ylabel('Normalized Weight of Each DOF')
subplot(1,2,2)
bar(LeftPC1,'DisplayName','LeftPC1')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 1                                                           .'))

figure()
subplot(1,2,1)
bar(RightPC2,'DisplayName','RightPC2')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Impaired Arm')
ylabel('Normalized Weight of Each DOF')
subplot(1,2,2)
bar(LeftPC2,'DisplayName','LeftPC2')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 2                                                           .'))

figure()
subplot(1,2,1)
bar(RightPC3,'DisplayName','RightPC3')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Impaired Arm')
ylabel('Normalized Weight of Each DOF')
subplot(1,2,2)
bar(LeftPC3,'DisplayName','LeftPC3')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 3                                                           .'))

figure()
subplot(1,2,1)
bar(RightPC4,'DisplayName','RightPC4')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Impaired Arm')
ylabel('Normalized Weight of Each DOF')
subplot(1,2,2)
bar(LeftPC4,'DisplayName','LeftPC4')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 4                                                           .'))

figure()
subplot(1,2,1)
bar(RightPC5,'DisplayName','RightPC5')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Impaired Arm')
ylabel('Normalized Weight of Each DOF')
subplot(1,2,2)
bar(LeftPC5,'DisplayName','LeftPC5')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 5                                                           .'))

figure()
subplot(5,2,1)
bar(RightPC1,'DisplayName','RightPC1')
axis([0 11 -1 1])
subplot(5,2,2)
bar(LeftPC1,'DisplayName','LeftPC1')
axis([0 11 -1 1])
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 1                                                           .'))

subplot(5,2,3)
bar(RightPC2,'DisplayName','RightPC2')
axis([0 11 -1 1])
subplot(5,2,4)
bar(LeftPC2,'DisplayName','LeftPC2')
axis([0 11 -1 1])
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 2                                                           .'))


subplot(5,2,5)
bar(RightPC3,'DisplayName','RightPC3')
axis([0 11 -1 1])
ylabel('Normalized Weight of Each DOF')
subplot(5,2,6)
bar(LeftPC3,'DisplayName','LeftPC3')
axis([0 11 -1 1])
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 3                                                           .'))

subplot(5,2,7)
bar(RightPC4,'DisplayName','RightPC4')
axis([0 11 -1 1])
subplot(5,2,8)
bar(LeftPC4,'DisplayName','LeftPC4')
axis([0 11 -1 1])
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 4                                                           .'))

subplot(5,2,9)
bar(RightPC5,'DisplayName','RightPC5')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Impaired Arm')
subplot(5,2,10)
bar(LeftPC5,'DisplayName','LeftPC5')
axis([0 11 -1 1])
xlabel('Kinematic DOFs of Unimpaired Arm')
title(strcat('SS', num2str(subjectID),' Kinematic Synergy Vector 5                                                           .'))

