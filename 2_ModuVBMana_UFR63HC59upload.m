%% 1. select the T1 image for GHR63 HC59, add computed TIV data to the Behavioral data;% UFR/HC
dataA = importdata('TIV.txt');
dataAta = array2table(dataA(:, [1:4]),'VariableNames', {'TIV', 'GM',...
    'WM', 'CSF'});
HCdataMatchImageFDTIVspm1 = [HCdataMatchImageFDTIVspm(:, [1:10]), dataAta,...
    HCdataMatchImageFDTIVspm(:, [15: end])];

%% compute the volume of regional putamen
% step1: obtain the ROI seed(subarea of striatum in AAL3);
path = ['~/3ndComputVBMcat12modulated'];
mask1 = zeros(113,137,113);
V1 = spm_vol([path, '/slice_AAL3v1.nii']);% 1.5*1.5*1.5;
T1 = spm_read_vols(V1);
for i = 77:78
    maskX = mask1+ (T1==i);
    sliceN = ['slice_AAL3striatum', num2str(i), 'mask.nii']
    V1.fname =  [path, filesep, sliceN];
    mask2 = spm_write_vol(V1, maskX);
end
% step2:extract volume of the subarea of striatum 
workpath = ['~/3ndComputVBMcat12modulated'];
groupP = [workpath, '/GHR63T1/mri'];%%%%%%%GHR63/HC59;
subfolder = dir([groupP, '/s6mwp1sub-*_T1w.nii']);

volume = zeros(113,137,113); 
roiVol = zeros(length(subfolder), 3);

modicell = {'slice_AAL3striatum77mask.nii', 'slice_AAL3striatum78mask.nii'};
clear roiVol;
for i = 1:length(subfolder)
    splitcell = split(subfolder(i).name, '_');
    split1cell = split(splitcell{1}, '-');
    roiVol(i, 1) = str2num(split1cell{2});
    Vsub = spm_vol([groupP, '/', subfolder(i).name]);
    volume1 = spm_read_vols(Vsub);
    
    for roii=1:length(modicell)
        clear plus;
        Vroi = spm_vol([workpath, '/', modicell{roii}]);
        volume = spm_read_vols(Vroi);
        plus = volume1.*volume;
        total = sum(plus(:));
        roiVol(i, roii+1) = total* 1.5^3/10^3;% mm3 --> ml
    end
end
roiVoltable = array2table(roiVol, 'VariableNames',...
    {'ID', 'put77volml', 'put78colml'});
    
%% %%%%subgroup the GHR by striatum(including putamen...) volume-median from HC;
%% compute the mean(SD) striatum(including putamen...) volume of HC59;
colu = 7;
meansd(1,1) = nanmean(HC59subareaStriatumVOL1{:, colu});
meansd(1,2) = nanstd(HC59subareaStriatumVOL1{:, colu});
meansd(1,3) = median(HC59subareaStriatumVOL1{:, colu});
% GHR63
meansd(2,1) = nanmean(GHR63subareaStriatumVOL1{:, colu});
meansd(2,2) = nanstd(GHR63subareaStriatumVOL1{:, colu});
meansd(2,3) = median(GHR63subareaStriatumVOL1{:, colu});

%% %%%%%%%%%%%%%%%1. separate the GHR by cutoff defined as the median of volume in HC;
clear num;
cutof =0.6681; % vol median;
col = 7; % col=2 ---L_caudate;

num(1, 1) = sum(GHR63subareaStriatumVOL1{:, col}>cutof)
num(2, 1) = sum(GHR63subareaStriatumVOL1{:, col}<=cutof)

num(1, 2) = sum(HC59subareaStriatumVOL1{:, col}>cutof)
num(2, 2) = sum(HC59subareaStriatumVOL1{:, col}<=cutof)

%% 1) subgroup the GHR by putamen  volume-median of HCs;
% load GHR and HC volume table;
load('GHR63subareaStriatumVOL1.mat');
load('HC59subareaStriatumVOL1.mat');

% (1) putamen-median as cutoff, separately;
cutof = 0.6681;% vol median mean;
col = 3;% col=2 ---L_caudate;

GHRmaxindex78 = find(GHR63subareaStriatumVOL1{:, col}>cutof);
save GHRmaxHCmedian_index78.mat GHRmaxindex78; 
GHRminindex78 = find(GHR63subareaStriatumVOL1{:, col}<=cutof);
save GHRminHCmedian_index78.mat GHRminindex78; 

load('GHRmaxHCmedian_index78.mat');
load('GHRminHCmedian_index78.mat');

%% mean(SD) of subarea subgroup
colu = [5, 9, 11:14];
row1 = GHRmaxindex78;
row2 = GHRminindex78;
rowHC = [1:59];
clear meansd;
for i = 1:length(colu)
    % group;
    meansd(i,1) = nanmean(GHRdataMatchImageFDTIV1{row1, colu(i)});
    meansd(i,2) = nanstd(GHRdataMatchImageFDTIV1{row1, colu(i)});
    meansd(i,3) = nanmean(GHRdataMatchImageFDTIV1{row2, colu(i)});
    meansd(i,4) = nanstd(GHRdataMatchImageFDTIV1{row2, colu(i)});
    meansd(i,5) = nanmean(HCdataMatchImageFDTIVspm1{rowHC, colu(i)});%%%%HC
    meansd(i,6) = nanstd(HCdataMatchImageFDTIVspm1{rowHC, colu(i)});
end

%% anova WITH covariates--2023(3 group + hoc-post);
%  behavioral variable, symptom, cognitive,,,
rowHC = [1:59];
colu1 =[12:14];
row1 = GHRmaxindex78;%
row2 = GHRminindex78;
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9, 11]};...
    GHRdataMatchImageFDTIV1{row2, [5, 4, 9, 11]};...
    HCdataMatchImageFDTIVspm1{rowHC, [5, 4, 9, 11]}];
clear fp;
clear tbl11;
for i = 1:length(colu1)
    data2comp = [GHRdataMatchImageFDTIV1{row1, colu1(i)};...
        GHRdataMatchImageFDTIV1{row2, colu1(i)};...
        HCdataMatchImageFDTIVspm1{rowHC, colu1(i)}];
    grouplabel = [ones(length(row1), 1);2*ones(length(row2), 1);...
        3*ones(length(rowHC), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu,tiv
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
    
    [results,~,~,gnames] = multcompare(stats,"Dimension",[1]);
    tbl11{i, 1} = array2table(results,'VariableNames',...
    {'GroupA','GroupB' ,'LowerLimi','A_B','UpperLimit','Pvalue'});
end

%% %%%%%% GHR63-HC59 VBM Analysis;%%%%%%%%%%%%%%%%%%%%%
%% %%% volume analysis between each GHR63 group and HC59 %%%
%% 1. to extract the cluster volume in the group of GHR63 and HC59, separately;
%%%%%%%%%%%%%%% GHR63 HC59 VBM 
WKD = ['~/3ndComputVBMcat12modulated/ghr63HC59volume2ndanalysis'];
cellname = {'GHR63T1', 'HC59T1'};
clear cellclumask;
for cl = 1: 10
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1-10
clear GHRHCclusterVolnanmean;
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);% image file;
    for i = 1: length(cellclumask) % read mask;
        maskvol = spm_vol([WKD, '/twosamTtestagesexeduTIV/ClusterSaveGRFvoxel0001cluster005/',...
            cellclumask{i}]);% need revise;
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            disp(['Processing the ', num2str(subi), '-th sub of ', cellname{j}, '!']);
            splitcell = split(subfiles(subi).name, '-');
            split1cell = split(splitcell{2}, '_');
            signal(subi, 1) = str2num(split1cell{1});% save the sub ID;
            
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image

            con1_data = con1(maskAAL>0); %  this is a vector
            mean = nanmean(con1_data(:));
            signal(subi, i+1) = mean* 1.5^3/10^3;% mm3 
        end
    end
    GHR63HC59clusterVolmean{j, 1} = signal; 
end

% save
save  GHR63HC59cluster19VolmeanXin.mat GHR63HC59clusterVolmean;

%% 2.correcction between cluster volume and symptoms scores(SIPS/GAF); 
% using age,sex,edu,TIV as covariates;
% rowINDEX  = [1: 9, 11:63]; % Remove GHR10181;
rowINDEX  = [1:63];
% col = [27,51];% SIPS/GAF;
% col = [ 55];% LES;
col = [64, 65, 68, 72, 74, 76];
ROIco = [2: 20];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHR63HC59clusterVolmean{1, 1}(rowINDEX, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

%% %%%%%%%%%%%%% GHR subgroup and HC59 %%%%%%%%%%%%%%%%%%%%%
%% extract the subgroup age,sex,edu, TIV behavior data and image data;
%% two T using age.sex.edu,TIV as COVs compare the group of GHRmax-HC and GHRmin-HC;
% using GRF corrected the T statistical results and using the saved cluster
% mask of passing the GRF correction 
%% 1. to extract the cluster volume in the group of GHRmax and GHRmin, separately;
%%%%%%%%%%%%%%%  VBM 
WKD = ['~/3ndComputVBMcat12modulated/subgroup_by_HC7778vol_analysis'];
cellname = {'78RGHRmaxsubG', 'HC59T1'};
clear cellclumask;
for cl = 1: 12
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1-12
clear GHRmaxHCclusterVolml;
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);% image file;
    for i = 1: length(cellclumask) % read mask;
        maskvol = spm_vol([WKD, '/78GHRmax_HC59twoTagesexeduTIV/ClusterSaveGRFvoxel0001cluster005/',...
            cellclumask{i}]);% need revise;
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            disp(['Processing the ', num2str(subi), '-th sub of ', cellname{j}, '!']);
            splitcell = split(subfiles(subi).name, '-');
            split1cell = split(splitcell{2}, '_');
            signal(subi, 1) = str2num(split1cell{1});% save the sub ID;
            
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image
            
            con1_data = con1(maskAAL>0); %  this is a vector
            mean = nanmean(con1_data(:));
            signal(subi, i+1) = mean* 1.5^3/10^3;% mm3 --> ml, %%%%%%%% mean volume;
        end
    end
    GHRmaxHCclusterVolml{j, 1} = signal; 
end

%% 2.correcction between cluster volume and symptoms scores(SIPS/GAF); 
% using age,sex,edu,TIV as covariates;
% MAX group;
rowINDEX  = GHRmaxindex78;
col = [27,51];% SIPS/GAF, sickN;
% col = [55];% LES;
ROIco = [2: 13];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRmaxHCclusterVolml{1, 1}(:, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

% MIN group;
% rowINDEX  = GHRminindex78([1:5, 7:end]);% remove 10181 for abnormal LES;
rowINDEX  = GHRminindex78;
col = [27, 51];% SIPS/GAF, sickN;
% col = [55];% LES;
ROIco = [2: 4];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRminHCclusterVolml{1, 1}(:, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

% HC59 -MINHC cluster group;
rowINDEX  = [1:11, 13: 59];
col = [26, 27];% LES;
ROIco = [2: 4];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = HCdataMatchImageFDTIVspm1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRminHCclusterVolml{2, 1}(rowINDEX, ROIco(i)),...
            HCdataMatchImageFDTIVspm1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

clear fdrQ;
fdrQ = mafdr(pvalue, 'BHFDR', true);

%% %%%%%% mediation analysis %%%%
%% subgroup cluster  GHR vol-LESN-symptom(SIPS);;
Rowin = GHRminindex78([1:5, 7:31, 33: end]);%%% index-subgroup;
Xa = GHRminHCclusterVolml{1, 1}(:,2);%  vol
M = zscore(Xa([1:5, 7:31, 33: end]));
X = zscore(GHRdataMatchImageFDTIV1{Rowin, 55});% LESN
Y = zscore(GHRdataMatchImageFDTIV1{Rowin, 27});% sips total
COV = GHRdataMatchImageFDTIV1{Rowin, [5, 4, 9, 11]};% age, sex, edu, TIV;

[paths, stats] = mediation(X, Y, M, 'boottop', 'bootsamples', 3000,...
    'covs', COV, 'plots', 'verbose',...
     'names', {'LESN' 'SIPStotal' 'LcaudateVol'});

 
 %% add analysis of enrichment for GHRmin and GHRmax subgroup;
 %% % read the image data Cell{sub37*1}--GHRmin-voxel mediation
 WKD = ['~/3ndComputVBMcat12modulated'];
maskvol = spm_vol([WKD, '/slice_AAL3v1nocerelluBin.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask

clear niifiles;
loadpath = [WKD, filesep, '/subgroup_by_HC7778vol_analysis/78RGHRminsubG'];
niifiles = dir([loadpath, '/*.nii']);
for j = 1:length(niifiles)
    clear V;
    clear volume;
    V = spm_vol([loadpath, filesep, niifiles(j).name]);
    [con1, XYZ] = spm_read_vols(V);  %  spm read image
    
    con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
    volumes{j, 1} = con1_data; % Cell{sub37*load1234}--382723*1double
end

[m, n] = size(volumes); %37*1
Rowin = GHRminindex78;%%% index-subgroup;
X = GHRdataMatchImageFDTIV1{Rowin, 55};% LESN,
Y = GHRdataMatchImageFDTIV1{Rowin, 27};% sips total
COV = GHRdataMatchImageFDTIV1{Rowin, [5, 4, 9, 11]};% age, sex, edu, TIV;

for k = 1: length(volumes{1,1}) % for each location--382723
    display(['Computing the ', num2str(k), '-th location point!']);
    clear subidata;
    for row = 1: m
        subidata(row, 1) = volumes{row, 1}(k);
    end
    LOCdata{k, 1} = subidata;
    % fitlme model 1
    subidata(isnan(subidata) == 1)=[];
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;

    dataROICov = [subidata(:, 1), X, Y];% volume ,LESN, sipS;
    indNan = find(sum(isnan(dataROICov), 2)>0); % find the row including nan;
   indnoNan = setdiff([1: length(Rowin)],  indNan); % 32
    dataROICov1 = zscore(dataROICov(indnoNan', :));
    
    [paths, stats] = mediation(dataROICov1(:, 2), dataROICov1(:, 3), dataROICov1(:, 1), 'boottop', 'bootsamples', 3000,...
        'covs', COV(indnoNan', :), 'verbose',...
        'names', {'LESN' 'SIPStotal' 'LcaudateVol'});
    close all;
    table1 = [stats.mean; stats.ste; stats.z; stats.ci(:, :, 1); stats.ci(:, :, 2); stats.p];
    table2 = array2table(table1, 'VariableNames', {'a', 'b', 'c1', 'c', 'ab'},...
        'RowNames', {'Coef', 'STE', 'Z', 'CI lb', 'CI ub', 'p'});
    statsSave{k, 1} = table2; % save the results;
    
    Cz(k, 1) =stats.z(4); % statistical z of variable c;
    Cp(k, 1) =stats.p(4); % statistical p of variable c(total effect);
    ABz(k, 1) =stats.z(5); % statistical z of variable ab;
    ABp(k, 1) =stats.p(5); % statistical p of variable ab(indirect path);
end
clearvars -except statsSave Cz Cp ABz ABp;
%Cz
V.fname = 'voxelMediationGHR78min_Cz.nii';
Fcon = con1;
statvalue = Cz;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= statvalue;
spm_write_vol(V, Fcon);
%ABz
V.fname = 'voxelMediationGHR78min_ABz.nii';
Fcon1 = con1;
statvalue1 = ABz;
Fcon1(maskAAL<1)= 0;
Fcon1(maskAAL>0)= statvalue1;
spm_write_vol(V, Fcon1);

%%  step1:extract the mean z value in brain ROI of cerellum part of AAL3;
% in AAL3, [1:94, 121:170]; of which no label equals to 35,36,81,82;
 path = ['~/3ndComputVBMcat12modulated'];
mask1 = zeros(113,137,113);
Vatlas = spm_vol([path, '/slice_AAL3v1.nii']);% 1.5*1.5*1.5;
Tatlas = spm_read_vols(Vatlas); % AAL3 data;

Vmap = spm_vol([path, '/subgroup_by_HC7778vol_analysis/addAnavoxelMediationGHRmin/voxelMediationGHR78min_ABz.nii']);% 1.5*1.5*1.5;
Tmap = spm_read_vols(Vmap); % Tmap data;
label = [1:94, 121:170];
for i = 1:length(label)
    maskX = mask1+ (Tatlas == label(i));
    data = Tmap(maskX > 0);
    voxelMediationGHR78minAAL3mean(i, 1) = label(i);
    voxelMediationGHR78minAAL3mean(i, 2) = nanmean(data);
end
% save the aal3 ROI144 mean Tvalue(144*2) as csv

%% %%%%%%%%%%%%GHRmax: similiar voxel mediation with above GHR78min subgroup%%%
%% % read the image data Cell{sub26*1}--GHRmax-voxel mediation
