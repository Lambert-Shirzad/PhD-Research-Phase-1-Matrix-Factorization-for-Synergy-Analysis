% This script takes as input the XLSX OpenSim output (converted from .mot)
% with the FEATHERS project OpenSim model's .mot template for the arrangement of data 
% (time series of joint positions). The main functionality is to find the columns of interest
% from among ~170 joints and save them in a .mat file for further analysis. 
% SubjectIDs will have the ID number of the subjects to be processed
% NumberPeriods has how many 40s periods subjects have, for example, if all 
% subjects have data for 0-320s, input NumberPeriods= 320/40 = 8.
% 20151209 Written by Navid Shirzad


function Xlsx2Mat_BatchProcess(SubjectIDs, NumberPeriods)


    for subjectcounter = 1:size(SubjectIDs,2)

        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end

        for periodcounter = 1:NumberPeriods
            TimePeriod = strcat(num2str((periodcounter-1)*40), '-', num2str(periodcounter*40), 's');
            xlsxDataFileName = strcat('Y', SubjID, '_', TimePeriod);
            % Read in the data 
            CurrentDirectory = cd;
            xlsxFileString = strcat(CurrentDirectory, '\', xlsxDataFileName, '.xlsx');
            [NUMERIC, TEXT] = xlsread(xlsxFileString); %Numeric now has all the data in it

            % find the column numbers of the joints of interest
            TIME = find(strcmp(TEXT(7,:), 'time'));
            backR1 = find(strcmp(TEXT(7,:), 'back_tilt'));
            backR2 = find(strcmp(TEXT(7,:), 'back_list'));
            backR3 = find(strcmp(TEXT(7,:), 'back_rotation'));
            RshoulderE = find(strcmp(TEXT(7,:), 'elv_angle'));
            RshoulderR1 = find(strcmp(TEXT(7,:), 'shoulder_elv'));
            RshoulderR2 = find(strcmp(TEXT(7,:), 'shoulder1_r2'));
            RshoulderR3 = find(strcmp(TEXT(7,:), 'shoulder_rot'));
            RelbowR1 = find(strcmp(TEXT(7,:), 'elbow_flexion'));
            RelbowR2 = find(strcmp(TEXT(7,:), 'pro_sup'));
            RwristR1 = find(strcmp(TEXT(7,:), 'deviation'));
            RwristR2 = find(strcmp(TEXT(7,:), 'flexion'));
            LshoulderE = find(strcmp(TEXT(7,:), 'elv_angle_l'));
            LshoulderR1 = find(strcmp(TEXT(7,:), 'shoulder_elv_l'));
            LshoulderR2 = find(strcmp(TEXT(7,:), 'shoulder1_r2_l'));
            LshoulderR3 = find(strcmp(TEXT(7,:), 'shoulder_rot_l'));
            LelbowR1 = find(strcmp(TEXT(7,:), 'elbow_flexion_l'));
            LelbowR2 = find(strcmp(TEXT(7,:), 'pro_sup_l'));
            LwristR1 = find(strcmp(TEXT(7,:), 'deviation_l'));
            LwristR2 = find(strcmp(TEXT(7,:), 'flexion_l'));

            % arrange the time serries that will be saved
            NumericData = [NUMERIC(:,TIME)...
                           NUMERIC(:,backR1) NUMERIC(:,backR2) NUMERIC(:,backR3)...
                           NUMERIC(:,RshoulderE) NUMERIC(:,RshoulderR1) NUMERIC(:,RshoulderR2) NUMERIC(:,RshoulderR3)...
                           NUMERIC(:,RelbowR1) NUMERIC(:,RelbowR2) NUMERIC(:,RwristR1) NUMERIC(:,RwristR2)...
                           NUMERIC(:,LshoulderE) NUMERIC(:,LshoulderR1) NUMERIC(:,LshoulderR2) NUMERIC(:,LshoulderR3)...
                           NUMERIC(:,LelbowR1) NUMERIC(:,LelbowR2) NUMERIC(:,LwristR1) NUMERIC(:,LwristR2)];
            HEADERS = {'time', 'back_tilt', 'back_list', 'back_rotation', 'elv_angle', 'shoulder_elv', ...
                       'shoulder1_r2', 'shoulder_rot', 'elbow_flexion', 'pro_sup', 'deviation', 'flexion', ...
                       'elv_angle_l', 'shoulder_elv_l', 'shoulder1_r2_l', 'shoulder_rot_l', 'elbow_flexion_l', ...
                       'pro_sup_l', 'deviation_l', 'flexion_l'};

            % save data in mat file
            OutputNum = strcat('Y', SubjID, '_', TimePeriod, '.mat');
            % OutputHeader = strcat('Headers_', SubjID, '.mat');
            % save(OutputHeader, 'HEADERS');
            save(OutputNum, 'NumericData');
        end
    end
end