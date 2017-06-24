% This script takes as input the XLSX Kinect output with the FEATHERS project's
% template for the arrangement of data (time series of position data for
% different joing positions). The main functionality is to resample the
% joint position data at a constant frequency and then filter the high
% frequency noise in the data. After Filtering, the orientation of the data 
% is changed from Kinect orientation to OpenSim data (OpenSim.X=Kinect.(-Z),
% OpenSim.Y=Kinect.Y, OpenSim.Z=Kinect.X). The resampled, filtered, reoriented
% data is then written in a seperate xlsx file.

% 20150321 Written by Navid Shirzad


function ResampleData(ResampleFreq, Cutoff_freq, xlsxDataFileName, SubjID)

% Read in the data required for resampling, i.e. the time stamps and all the position-related data.
CurrentDirectory = cd;
xlsxFileString = strcat(CurrentDirectory, '\', xlsxDataFileName, '.xlsx');
StartTimeStamp = xlsread(xlsxFileString, 2, 'D2'); %D2 in sheet 4 has the time stamp for the start of the trial
NUMERIC = xlsread(xlsxFileString,1); %Numeric now has all the data in it
StartRow = find(NUMERIC(:,3) > StartTimeStamp, 1);
NUMERIC = NUMERIC(StartRow:end, 3:end); %columns 1 and 2 are NaN

% define the column number of the marker positions, update if needed
TIME = 1;
HipCenter = 51;
ShoulderCenter = 39;
ShoulderLeft = 31;
ElbowLeft = 23;
WristLeft = 5;
HandLeft = 15;
ShoulderRight = 35;
ElbowRight = 27;
WristRight = 10;
HandRight = 19;


% arrange the time serries that will be resampled
DataToResample = [NUMERIC(:,HipCenter:HipCenter+2) NUMERIC(:,ShoulderCenter:ShoulderCenter+2)...
    NUMERIC(:,ShoulderLeft:ShoulderLeft+2) NUMERIC(:,ElbowLeft:ElbowLeft+2)...
    NUMERIC(:,WristLeft:WristLeft+2) NUMERIC(:,HandLeft:HandLeft+2)...
    NUMERIC(:,ShoulderRight:ShoulderRight+2) NUMERIC(:,ElbowRight:ElbowRight+2)...
    NUMERIC(:,WristRight:WristRight+2) NUMERIC(:,HandRight:HandRight+2)];
    %three columns have xyz position data in the kinect orientation
Original_t = NUMERIC(1:end,TIME)/1000; % convert to seconds 

% Shift the timestamps so it starts from 0 seconds
Shifted_t = [0;Original_t(2:end,1)-Original_t(1,1)];

% create a timeseries object from the original data
Original_timeseries = timeseries(DataToResample,Shifted_t); 

% create a timeline based on the desired resampling frequency
Desired_timeseries = [0:1/ResampleFreq:Shifted_t(end,1)]';

% resample the originaltimeseries
Resampled_data = resample(Original_timeseries,Desired_timeseries);

% filter the data using a butterworth filter
Filtered_data_wo_time = ButterworthFilter(ResampleFreq, Cutoff_freq, Resampled_data);
%Filtered_data = [Desired_timeseries, Filtered_data_wo_time]; 

% change the orientation
% (OpenSim.X=Kinect.(-Z),OpenSim.Y=Kinect.Y, OpenSim.Z=Kinect.X)
numMarkers = 10;
Oriented_data_wo_time = Filtered_data_wo_time;
for i=1:numMarkers
    temp = Oriented_data_wo_time(:,i*3-2);
    Oriented_data_wo_time(:,i*3-2) = -Oriented_data_wo_time(:,i*3);
    Oriented_data_wo_time(:,i*3) = temp;
end
Final_data = [Desired_timeseries, Oriented_data_wo_time];

% write data in xlsx file
OutputName = strcat(SubjID, '_resampled&filtered&reoriented_UE_markers.xls');
xlswrite(OutputName, Final_data);
end