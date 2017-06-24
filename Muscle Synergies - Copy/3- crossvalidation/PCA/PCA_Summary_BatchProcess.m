% to extract NNMF performance summary data 
% Outcome will be several .mat files

% 20160330 Written by Navid Shirzad

function PCA_Summary_BatchProcess(SubjectIDs)
    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end
        
        load(strcat('PCA_Subj_', SubjID, '.mat' ));
        
        if SubjDomHand == 'R'
            
            ndimDom(subjectcounter,1) = mean(ndimR);
            CAFDom(subjectcounter,1) = mean(CAFRightTr);
            AvgReconsErrTrDom(subjectcounter,1) = mean(AvgReconsErrRightTr);
            AvgReconsErrValDom(subjectcounter,1) = mean(AvgReconsErrRightVal);
            RegCoeffTrDom(subjectcounter,1) = mean(RegCoeffRightTr);
            RegCoeffValDom(subjectcounter,1) = mean(RegCoeffRightVal);
            nCommonDom(subjectcounter,1) = mean(nCommonR);
            
            ndimNonDom(subjectcounter,1) = mean(ndimL);
            CAFNonDom(subjectcounter,1) = mean(CAFLeftTr);
            AvgReconsErrTrNonDom(subjectcounter,1) = mean(AvgReconsErrLeftTr);
            AvgReconsErrValNonDom(subjectcounter,1) = mean(AvgReconsErrLeftVal);
            RegCoeffTrNonDom(subjectcounter,1) = mean(RegCoeffLeftTr);
            RegCoeffValNonDom(subjectcounter,1) = mean(RegCoeffLeftVal);
            nCommonNonDom(subjectcounter,1) = mean(nCommonL);
                        
            ElapsedTime(subjectcounter,1) = timeElapsed;
            nCommonBothSides(subjectcounter,1) = mean(nCommon);
            
        else
            
            ndimDom(subjectcounter,1) = mean(ndimL);
            CAFDom(subjectcounter,1) = mean(CAFLeftTr);
            AvgReconsErrTrDom(subjectcounter,1) = mean(AvgReconsErrLeftTr);
            AvgReconsErrValDom(subjectcounter,1) = mean(AvgReconsErrLeftVal);
            RegCoeffTrDom(subjectcounter,1) = mean(RegCoeffLeftTr);
            RegCoeffValDom(subjectcounter,1) = mean(RegCoeffLeftVal);
            nCommonDom(subjectcounter,1) = mean(nCommonL);
            
            ndimNonDom(subjectcounter,1) = mean(ndimR);
            CAFNonDom(subjectcounter,1) = mean(CAFRightTr);
            AvgReconsErrTrNonDom(subjectcounter,1) = mean(AvgReconsErrRightTr);
            AvgReconsErrValNonDom(subjectcounter,1) = mean(AvgReconsErrRightVal);
            RegCoeffTrNonDom(subjectcounter,1) = mean(RegCoeffRightTr);
            RegCoeffValNonDom(subjectcounter,1) = mean(RegCoeffRightVal);
            nCommonNonDom(subjectcounter,1) = mean(nCommonR);
                        
            ElapsedTime(subjectcounter,1) = timeElapsed;
            nCommonBothSides(subjectcounter,1) = mean(nCommon);
            
        end
    end
              
        %% save the analysis outcome
         
                     
            save(strcat('PCA_Performance_Summary', '.mat' ), ...
            'ndimDom', 'CAFDom', 'AvgReconsErrTrDom', 'AvgReconsErrValDom', ...
            'RegCoeffTrDom', 'RegCoeffValDom', 'nCommonDom', ...
            'ndimNonDom', 'CAFNonDom', 'AvgReconsErrTrNonDom', 'AvgReconsErrValNonDom', ...
            'RegCoeffTrNonDom', 'RegCoeffValNonDom', 'nCommonNonDom', ...
            'ElapsedTime', 'nCommonBothSides');
        
end