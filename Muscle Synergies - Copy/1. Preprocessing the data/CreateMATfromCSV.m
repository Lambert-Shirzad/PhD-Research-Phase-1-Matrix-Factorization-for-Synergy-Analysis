% This script takes in raw EMG data from 16 channels (8 channels from each
% side of the body) in CSV format, and saves the data in .mat format. 

%change lines 19,29,30 to go between Weighted and Lucky Pirate versions of
%the data.

% 20160321 Written by Navid Shirzad


function CreateMATfromCSV(SubjectIDs)
    
    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end
        
        xlsxDataFileName = strcat('RAW_Weight_Subj', SubjID, '.csv');
        CurrentDirectory = cd;
        xlsxFileString = strcat(CurrentDirectory, '\', xlsxDataFileName);
        [NUMERIC, TEXT] = xlsread(xlsxFileString); %Numeric now has all the data in it
        
        %time = NUMERIC(1:610001,1);
        rightside = NUMERIC(:,1:9);
        leftside = [NUMERIC(:,1) , NUMERIC(:,10:17)];

        % save data in mat file
        OutputRight = strcat('RAW_Weight_Subj_', SubjID, '_Right.mat');
        OutputLeft = strcat('RAW_Weight_Subj_', SubjID, '_Left.mat');
        %HEADERS = {'time','Delt_Ant','Delt_Med','Delt_Post','Biceps','Triceps_Long','Triceps_Lat','Brachi','Pect-Maj'};
        %save('HEADERS.mat', 'HEADERS');
        save(OutputRight, 'rightside');
        save(OutputLeft, 'leftside');
        
    end
    
