% This script combines k .mat files that have joint displacement data in 
% them to create the full set data for a subject. 
% Use this to batch process the data

% 20151209 Written by Navid Shirzad


function CreateFullSet(SubjectIDs, NumberPeriods)

Duration = 40;

load('SS3_0-40s.mat')
FullSet = NumericData(1,:);

%there are k=10 folds

for k = 1:7
    %load the mat files one by one
    load(strcat('SS3_', num2str((k-1)*Duration), '-', num2str((k)*Duration), 's.mat' ));
    FullSet = [FullSet; NumericData(2:end,:)];
end


% % % % % % % %use the following lines if the time period of the last mat file is not
% % % % % % % %exactly 40 seconds. make sure the subject number is updated.
% % % % % % % 
% % % % % % % load('SS3_280-315.1s.mat')
% % % % % % % FullSet = [FullSet; NumericData(2:end,:)];

save('FullSet_SS3.mat', 'FullSet')


