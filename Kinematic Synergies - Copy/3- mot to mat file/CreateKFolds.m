% This script combines k .mat files that have joint displacement data in 
% them to create k-folds of training and evaluation .mat files. Evaluation 
% files are exactly the same old .mat files, and training files have all 
% the data from k-1 .mat files in them.

% 20150330 Written by Navid Shirzad

%each fold has Duration = 40 seconds of data
Duration = 40;

load('SS4_0-40s.mat')
FirstLine = NumericData(1,:);
%there are k=10 folds
for kprime=1:10
    TrainingNumeric = FirstLine;
    for k = 1:10
        %load the mat files one by one
        if k ~= kprime 
            load(strcat('SS4_', num2str((k-1)*Duration), '-', num2str((k)*Duration), 's.mat' ));
            TrainingNumeric = [TrainingNumeric; NumericData(2:end,:)];
        end
    end
    save(strcat('SS4_Training_Exclude_', num2str((kprime-1)*Duration), '-', num2str((kprime)*Duration), 's.mat' ), 'TrainingNumeric')
    clear TrainingNumeric
end

