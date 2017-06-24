% to extract mix all performance summary .mat files of different methods
% into one xlsx file

% 20160330 Written by Navid Shirzad

function mat2xlsx_SummaryGenerator()
    
    load('PCA_Performance_Summary.mat');
    sizePCA = size(ndimDom,1);%how many data points are there?
    METHODS(1:sizePCA*2,1) = 1; %1 means PCA, multiplied by2 b/c there are two hands.
    HandDom(1:sizePCA,1) = 1; %dominant hand
    HandDom(sizePCA+1:sizePCA*2,1) = 2; %non-dominant hand
    ndimAll(1:sizePCA,1) = ndimDom;
    ndimAll(sizePCA+1:sizePCA*2,1) = ndimNonDom;
    AvgReconsErrTrAll(1:sizePCA,1) = AvgReconsErrTrDom;
    AvgReconsErrTrAll(sizePCA+1:sizePCA*2,1) = AvgReconsErrTrNonDom;
    AvgReconsErrValAll(1:sizePCA,1) = AvgReconsErrValDom;
    AvgReconsErrValAll(sizePCA+1:sizePCA*2,1) = AvgReconsErrValNonDom;
    
    for i = 1:sizePCA
        if RegCoeffTrDom(i,1) > 1 
            RegCoeffTrDom(i,1) = 1/RegCoeffTrDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffTrNonDom(i,1) > 1 
            RegCoeffTrNonDom(i,1) = 1/RegCoeffTrNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValDom(i,1) > 1 
            RegCoeffValDom(i,1) = 1/RegCoeffValDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValNonDom(i,1) > 1 
            RegCoeffValNonDom(i,1) = 1/RegCoeffValNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
    end
    
    RegCoeffTrAll(1:sizePCA,1) = RegCoeffTrDom;
    RegCoeffTrAll(sizePCA+1:sizePCA*2,1) = RegCoeffTrNonDom;
    RegCoeffValAll(1:sizePCA,1) = RegCoeffValDom;
    RegCoeffValAll(sizePCA+1:sizePCA*2,1) = RegCoeffValNonDom;
    nCommonAll(1:sizePCA,1) = nCommonDom;
    nCommonAll(sizePCA+1:sizePCA*2,1) = nCommonNonDom;
    
    nCommonBothSidesAll(1:sizePCA,1) = 1; %1 means it is PCA
    nCommonBothSidesAll(1:sizePCA,2) = nCommonBothSides;
    ElapsedTimeAll(1:sizePCA,1) = ElapsedTime;
    
    load('NNMF_Performance_Summary.mat');
    sizeNNMF = size(ndimDom,1);%how many data points are there?
    METHODS(sizePCA*2+1:sizePCA*2+sizeNNMF*2,1) = 2; %2 means NNMF, multiplied by2 b/c there are two hands.
    HandDom(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = 1; %dominant hand
    HandDom(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = 2; %non-dominant hand
    ndimAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = ndimDom;
    ndimAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = ndimNonDom;
    AvgReconsErrTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = AvgReconsErrTrDom;
    AvgReconsErrTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = AvgReconsErrTrNonDom;
    AvgReconsErrValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = AvgReconsErrValDom;
    AvgReconsErrValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = AvgReconsErrValNonDom;
    
    for i = 1:sizeNNMF
        if RegCoeffTrDom(i,1) > 1 
            RegCoeffTrDom(i,1) = 1/RegCoeffTrDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffTrNonDom(i,1) > 1 
            RegCoeffTrNonDom(i,1) = 1/RegCoeffTrNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValDom(i,1) > 1 
            RegCoeffValDom(i,1) = 1/RegCoeffValDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValNonDom(i,1) > 1 
            RegCoeffValNonDom(i,1) = 1/RegCoeffValNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
    end
    
    RegCoeffTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = RegCoeffTrDom;
    RegCoeffTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = RegCoeffTrNonDom;
    RegCoeffValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = RegCoeffValDom;
    RegCoeffValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = RegCoeffValNonDom;
    nCommonAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1) = nCommonDom;
    nCommonAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1) = nCommonNonDom;    
    
    nCommonBothSidesAll(sizePCA+1:sizePCA+sizeNNMF,1) = 2; %2 means it is NNMF
    nCommonBothSidesAll(sizePCA+1:sizePCA+sizeNNMF,2) = nCommonBothSides;
    ElapsedTimeAll(sizePCA+1:sizePCA+sizeNNMF,1) = ElapsedTime;
    
    load('ICA_Performance_Summary.mat');
    sizeICA = size(ndimDom,1);%how many data points are there?
    METHODS(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = 3; %3 means ICA, multiplied by2 b/c there are two hands.
    HandDom(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = 1; %dominant hand
    HandDom(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = 2; %non-dominant hand
    ndimAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = ndimDom;
    ndimAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = ndimNonDom;
    AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = AvgReconsErrTrDom;
    AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = AvgReconsErrTrNonDom;
    AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = AvgReconsErrValDom;
    AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = AvgReconsErrValNonDom;
    
    for i = 1:sizeICA
        if RegCoeffTrDom(i,1) > 1 
            RegCoeffTrDom(i,1) = 1/RegCoeffTrDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffTrNonDom(i,1) > 1 
            RegCoeffTrNonDom(i,1) = 1/RegCoeffTrNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValDom(i,1) > 1 
            RegCoeffValDom(i,1) = 1/RegCoeffValDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
        if RegCoeffValNonDom(i,1) > 1 
            RegCoeffValNonDom(i,1) = 1/RegCoeffValNonDom(i,1);%invert the slope, let's keep all the slopes below 1. this should be a measure of centrality, not how far below or above of 1 it is. 
        end
    end
    
    RegCoeffTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = RegCoeffTrDom;
    RegCoeffTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = RegCoeffTrNonDom;
    RegCoeffValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = RegCoeffValDom;
    RegCoeffValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = RegCoeffValNonDom;
    nCommonAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1) = nCommonDom;
    nCommonAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1) = nCommonNonDom;
    
    nCommonBothSidesAll(sizePCA+sizeNNMF+1:sizePCA+sizeNNMF+sizeICA,1) = 3; %3 means it is ICA
    nCommonBothSidesAll(sizePCA+sizeNNMF+1:sizePCA+sizeNNMF+sizeICA,2) = nCommonBothSides;
    ElapsedTimeAll(sizePCA+sizeNNMF+1:sizePCA+sizeNNMF+sizeICA,1) = ElapsedTime;
    
    %save the data
   
%     OutputName = 'Factorization_Performance_Metrics_Healthy_Kinematic_Synergies.xls';
%     xlswrite(OutputName, 'Columns_are:_PCA=1/NNMF=2/ICA=3 and 1 for dominant hand and 2 for nondominant hand', 1, 'A1');
%     xlswrite(OutputName, 'Sheets are: ndimAll AvgReconsErrTrAll AvgReconsErrValAll RegCoeffTrAll RegCoeffValAll nCommonAll', 1, 'A2');
%     xlswrite(OutputName, [METHODS, HandDom, ndimAll, AvgReconsErrTrAll, AvgReconsErrValAll, RegCoeffTrAll, RegCoeffValAll, nCommonAll], 2 );
%     xlswrite(OutputName, [nCommonBothSidesAll, ElapsedTimeAll], 3);
    
    %plot the data
    
% %     figure()
% %     boxplot([ndimAll(1:sizePCA,1), ndimAll(sizePCA+1:sizePCA*2,1), ...
% %         ndimAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), ndimAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         ndimAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), ndimAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0 10]);
% %     ylabel('# Identified Synergy Vectors');
% %     
% %     figure()
% %     boxplot([AvgReconsErrTrAll(1:sizePCA,1), AvgReconsErrTrAll(sizePCA+1:sizePCA*2,1), ...
% %         AvgReconsErrTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), AvgReconsErrTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0 20]);
% %     ylabel('Average Reconstruction Error of Training Datasets');
% %     
% %     figure()
% %     boxplot([AvgReconsErrValAll(1:sizePCA,1), AvgReconsErrValAll(sizePCA+1:sizePCA*2,1), ...
% %         AvgReconsErrValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), AvgReconsErrValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0 20]);
% %     ylabel('Average reconstruction Error of Validation Datasets');
% %     
% %     figure()
% %     boxplot([RegCoeffTrAll(1:sizePCA,1), RegCoeffTrAll(sizePCA+1:sizePCA*2,1), ...
% %         RegCoeffTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), RegCoeffTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         RegCoeffTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), RegCoeffTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0.8 1]);
% %     ylabel('Correlation Coefficient of Original & Reconstructed Training Datasets');
% %     
% %     figure()
% %     boxplot([RegCoeffValAll(1:sizePCA,1), RegCoeffValAll(sizePCA+1:sizePCA*2,1), ...
% %         RegCoeffValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), RegCoeffValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         RegCoeffValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), RegCoeffValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0.8 1]);
% %     ylabel('Correlation Coefficient of Original & Reconstructed Validation Datasets');
% %             
% %     figure()
% %     boxplot(100*[nCommonAll(1:sizePCA,1), nCommonAll(sizePCA+1:sizePCA*2,1), ...
% %         nCommonAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), nCommonAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
% %         nCommonAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), nCommonAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
% %     axis([0.5 6.5 0 120]);
% %     ylabel('% Common Synergy Vectors btw Training & Validation Datasets');
% %     
% %     figure()
% %     boxplot(100*[nCommonBothSidesAll(1:sizePCA,2), ...
% %         nCommonBothSidesAll(sizePCA+1:sizePCA+sizeNNMF,2), ...
% %         nCommonBothSidesAll(sizePCA+sizeNNMF+1:sizePCA+sizeNNMF+sizeICA,2)], ...
% %         'labels', {'PCA','NNMF','ICA'})
% %     axis([0.5 3.5 0 150]);
% %     ylabel('% Common Synergy Vectors btw Right & Left Limb');
% %     
% %     figure()
% %     boxplot([ElapsedTimeAll(1:sizePCA,1), ...
% %         ElapsedTimeAll(sizePCA+1:sizePCA+sizeNNMF,1), ...
% %         ElapsedTimeAll(sizePCA+sizeNNMF+1:sizePCA+sizeNNMF+sizeICA,1)], ...
% %         'labels', {'PCA','NNMF','ICA'})
% %     axis([0.5 3.5 0 35]);
% %     ylabel('Training Time (s)');
    
    figure()
    subplot(3,2,1)
    boxplot([ndimAll(1:sizePCA,1), ndimAll(sizePCA+1:sizePCA*2,1), ...
        ndimAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), ndimAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        ndimAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), ndimAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0 10]);
    ylabel('# Identified Synergy Vectors');
    
    subplot(3,2,3)
    boxplot([AvgReconsErrTrAll(1:sizePCA,1), AvgReconsErrTrAll(sizePCA+1:sizePCA*2,1), ...
        AvgReconsErrTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), AvgReconsErrTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), AvgReconsErrTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0 20]);
    ylabel({'Average Reconstruction Error', 'of Training Datasets'});
    
    subplot(3,2,5)
    boxplot([AvgReconsErrValAll(1:sizePCA,1), AvgReconsErrValAll(sizePCA+1:sizePCA*2,1), ...
        AvgReconsErrValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), AvgReconsErrValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), AvgReconsErrValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0 20]);
    ylabel({'Average reconstruction Error', 'of Validation Datasets'});
    
    subplot(3,2,2)
    boxplot([RegCoeffTrAll(1:sizePCA,1), RegCoeffTrAll(sizePCA+1:sizePCA*2,1), ...
        RegCoeffTrAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), RegCoeffTrAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        RegCoeffTrAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), RegCoeffTrAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0.8 1]);
    ylabel({'Correlation Coefficient of Original', '& Reconstructed Training Datasets'});
    
    subplot(3,2,4)
    boxplot([RegCoeffValAll(1:sizePCA,1), RegCoeffValAll(sizePCA+1:sizePCA*2,1), ...
        RegCoeffValAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), RegCoeffValAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        RegCoeffValAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), RegCoeffValAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0.8 1]);
    ylabel({'Correlation Coefficient of Original', '& Reconstructed Validation Datasets'});
            
    subplot(3,2,6)
    boxplot(100*[nCommonAll(1:sizePCA,1), nCommonAll(sizePCA+1:sizePCA*2,1), ...
        nCommonAll(sizePCA*2+1:sizePCA*2+sizeNNMF,1), nCommonAll(sizePCA*2+sizeNNMF+1:sizePCA*2+sizeNNMF*2,1), ...
        nCommonAll(sizePCA*2+sizeNNMF*2+1:sizePCA*2+sizeNNMF*2+sizeICA,1), nCommonAll(sizePCA*2+sizeNNMF*2+sizeICA+1:sizePCA*2+sizeNNMF*2+sizeICA*2,1)])
    axis([0.5 6.5 0 120]);
    ylabel({'% Common Synergy Vectors btw', 'Training & Validation Datasets'});
    
end