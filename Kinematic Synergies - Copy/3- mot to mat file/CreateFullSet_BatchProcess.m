% This script combines k .mat files that have joint displacement data in 
% them to create the full set data for a subject. 
% Use this to batch process the data

% 20151209 Written by Navid Shirzad


function CreateFullSet_BatchProcess(SubjectIDs, NumberPeriods)
    FullSet = zeros(1200*NumberPeriods,20); %1200 = 30Hz * 40s
    for subjectcounter = 1:size(SubjectIDs,2)

            if SubjectIDs(subjectcounter) < 10
                SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
            else
                SubjID = num2str(SubjectIDs(subjectcounter));
            end

            for periodcounter = 1:NumberPeriods
                TimePeriod = strcat(num2str((periodcounter-1)*40), '-', num2str(periodcounter*40), 's');
                MatFileName = strcat('Y', SubjID, '_', TimePeriod, '.mat');
                % Read in the data 
                load(MatFileName)
                FullSet((periodcounter-1)*1200+1:periodcounter*1200,:) = NumericData(2:end,:); %loosing data of t=0s
            end
            
            save(strcat('FullSet_', 'Y', SubjID, '.mat'), 'FullSet')

    end
end


