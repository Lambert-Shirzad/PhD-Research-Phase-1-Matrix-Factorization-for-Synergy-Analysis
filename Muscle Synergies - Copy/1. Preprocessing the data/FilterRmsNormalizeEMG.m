% This script takes in raw EMG data from 8 channels (one side of the body)
% in .mat format and filters, calculates RMS, and normalizes the data and
% saves the data in .mat format. 

%change lines 21,24,138,139 to go between Weighted and Lucky Pirate versions of
%the data.
% 20160322 Written by Navid Shirzad


function FilterRmsNormalizeEMG(SubjectIDs)
    
    for subjectcounter = 1:size(SubjectIDs,2)

        %% load the data
        if SubjectIDs(subjectcounter) < 10
            SubjID = strcat('0', num2str(SubjectIDs(subjectcounter)));
        else
            SubjID = num2str(SubjectIDs(subjectcounter));
        end
        
        DataFileNameRight = strcat('RAW_Subj_', SubjID, '_Right.mat');
        CurrentDirectory = cd;
        load(strcat(CurrentDirectory, '\', DataFileNameRight));    %loads a variable called rightside, EMGs of 8 muscles on the right side  
        DataFileNameLeft = strcat('RAW_Subj_', SubjID, '_Left.mat');
        load(strcat(CurrentDirectory, '\', DataFileNameLeft));    %loads a variable called leftside, EMGs of 8 muscles on the left side  

        %% design the filter and filter the data        
        FilterObject = designfilt('bandpassiir', ...       % Response type
        'StopbandFrequency1',12, ...    % Frequency constraints
        'PassbandFrequency1',20, ...
        'PassbandFrequency2',400, ...
        'StopbandFrequency2',450, ...
        'StopbandAttenuation1',40, ...   % Magnitude constraints
        'PassbandRipple',3, ...
        'StopbandAttenuation2',40, ...
        'DesignMethod','butter', ...      % Design method
        'SampleRate',2000);            % Sample rate
               
        
        FilteredRightSide = filter(FilterObject, rightside(:,2:9));    
        FilteredLeftSide = filter(FilterObject, leftside(:,2:9));   
        
        %% Calculate a moving window rms
        
        windowlength = 40; %data is 2kHz, so moving rms will be 100Hz (40 and 20)
        overlap = 20;
        delta = windowlength - overlap;
        
        indices = 1:delta:length(FilteredRightSide);
        if length(FilteredRightSide) - indices(end) + 1 < windowlength
            indices = indices(1:find(indices+windowlength-1 <= length(FilteredRightSide), 1, 'last'));
        end

        ProcessedRightSide = zeros(length(indices),size(FilteredRightSide,2)); %data has 8 columns
        ProcessedLeftSide = zeros(length(indices),size(FilteredLeftSide,2));
        % Square the samples
        FilteredRightSide = FilteredRightSide.^2;
        FilteredLeftSide = FilteredLeftSide.^2;

        index = 0;
        for i = indices
            index = index+1;
            % Average and take the square root of each window
            ProcessedRightSide(index,:) = sqrt(mean(FilteredRightSide(i:i+windowlength-1,:)));
            ProcessedLeftSide(index,:) = sqrt(mean(FilteredLeftSide(i:i+windowlength-1,:)));
        end
        
        %% Normalize each channel (column) to its maximum and save as a percentage
        
        maxRight = max(ProcessedRightSide);
        maxLeft = max(ProcessedLeftSide);
        
        for i=1:size(FilteredRightSide,2) %assuming right and left side have the same size
            ProcessedRightSide(:,i) = 100*ProcessedRightSide(:,i)/maxRight(1,i);
            ProcessedLeftSide(:,i) = 100*ProcessedLeftSide(:,i)/maxLeft(1,i);
            
        end
            
        %% Visualize the data
        
%         figure()  
%         plot(rightside(:,1), rightside(:,2));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,1));
%         title('Delt Ant')
% %         figure()  
% %         plot(indices/2000, ProcessedRightSide2(:,1));
% %         title('Delt Ant')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,3));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,2));
%         title('Delt Med')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,4));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,3));
%         title('Delt Post')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,5));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,4));
%         title('Biceps')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,6));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,5));
%         title('triceps Long')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,7));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,6));
%         title('triceps Lat')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,8));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,7));
%         title('Brachi')
%         
%         figure()  
%         plot(rightside(:,1), rightside(:,9));
%         hold on
%         plot(indices/2000, ProcessedRightSide(:,8));
%         title('Pect Maj')
        
        %% Save data in mat file
        
        %add the time stamps
        ProcessedRightSide(:,1:9) = [(indices/2000)'-0.0005, ProcessedRightSide(:,1:8)];
        ProcessedLeftSide(:,1:9) = [(indices/2000)'-0.0005, ProcessedLeftSide(:,1:8)];
        
        OutputRight = strcat('Processed_Subj_', SubjID, '_Right.mat');
        OutputLeft = strcat('Processed_Subj_', SubjID, '_Left.mat');
        %HEADERS = {'time','Delt_Ant','Delt_Med','Delt_Post','Biceps','Triceps_Long','Triceps_Lat','Brachi','Pect-Maj'};
        %save('HEADERS.mat', 'HEADERS');
        save(OutputRight, 'ProcessedRightSide');
        save(OutputLeft, 'ProcessedLeftSide');
        
    end
    
