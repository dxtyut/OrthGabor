clc;clear;close all;

%% Parameters
paraPatchSize = 17;
paraBlkSize = [8, 8];

FiltersName = {'OrthGabor','Gabor'};
GaborFilters =  getGabor(paraPatchSize,paraPatchSize);
figure;

%% Load Gallery
load('EYB_gallery.mat');
TrnData_ImgCell = mat2imgcell(double(TrnData),ImgSize(1),ImgSize(2),'gray'); 

%%% P is the mean patch and C is the ZCA whitening matrix
fprintf('====== Load or learn P and C  \n\n\n');
load(sprintf('P_C_PS%d.mat',paraPatchSize));
%%% [P,C]=learnZCA(TrnData_ImgCell,paraPatchSize);
    

for flagOrth = 1:2
%%% flagOrth = 1, use orth gabors
%%% flagOrth = 2, use gabors

    %% Filters
    fprintf('====== Get %s Filters \n',FiltersName{flagOrth});
    Filters = cell(1,5);
    for ii=1:5
        if flagOrth==1
            Filters{ii}=orth(GaborFilters((ii-1)*8+1:ii*8,:)');
        end
        if flagOrth==2
            Filters{ii}=GaborFilters((ii-1)*8+1:ii*8,:)';
        end
    end

    %% Training
    %%% get train features
    fprintf('====== Compute features for Gallery samples  \n');
    ftrain = [];
    for ii=1:length(Filters)
        ftrain = [ftrain;feaExtract(TrnData_ImgCell,Filters{ii},P,C, paraPatchSize,paraBlkSize)];
    end
    ftrain=normc(sqrt(ftrain));

    %% Testing 
    fprintf('====== Testing for image-544 of EYB_s5 \n');
    load('EYB_s5_image544.mat');
    TestData_ImgCell = mat2imgcell(TestData,ImgSize(1),ImgSize(2),'gray'); 

    %%% get test features
    ftest = [];
    for ii=1:length(Filters)
        ftest = [ftest;feaExtract(TestData_ImgCell,Filters{ii},P,C,paraPatchSize,paraBlkSize)];
    end
    ftest=normc(sqrt(ftest));

    %%  Classification  
    distance_trts = pdist2(ftest',ftrain','cosine');
    [~, label_estimate] = min(distance_trts,[],2);
    fprintf('====== True label 30, Predicted label %d for image-544 of EYB_s5 \n\n\n',label_estimate);
    
    %% Show Distances
    subplot(2,1,flagOrth);
    
    bar(distance_trts); xlabel('Gallery Image ID'); set(gca,'xtick',[1:1:38]); 
    if flagOrth==1
        ylim([0.75, 0.85]);
    else
        ylim([0.6, 0.65]);
    end    
    ylabel(['Distance between',sprintf('\n'),' this Probe and all Galleries']);
    box on;grid on;ax=gca;
    ax.FontName='Times New Roman';ax.GridLineStyle='-.';
    set(gca,'Linewidth',1.5);set(gca,'FontSize',15);
    
    pause(0.5);
end