% Analysis of abundance data from Yeast GFP Library from Dnervaud et al.
% 2013. PNAS. Raw data can be found at
% http://128.179.34.6/twiki/bin/view/CellImaging/WebHome.

clc
clear all


%% Load data

filename=load('MMS_high_info'); % An example with MMS_high_info.mat
StrName=fieldnames(filename);
StrName=StrName{1};

%% Find treatment timepoints

pre_treatment_time_ind = find(filename.(StrName).Time<0); % pre-treatment timepoints

treatment_time_ind = find(filename.(StrName).Time>0); % treatment timepoints

%% Disregard strains with abundance values of high variance

% Find variance and mean of abundance of every strain at pre-treatment timepoints

pre_treatment_BGminSubMean = filename.(StrName).BGminSubMean(:,pre_treatment_time_ind);

for i = 1:length(filename.(StrName).BGminSubMean)
    pre_treatment_variance(i,1) = var(pre_treatment_BGminSubMean(i,:));
    pre_treatment_mean(i,1) = mean(pre_treatment_BGminSubMean(i,:));
end

% Find strains that have a variance > 2-fold

strains_excluded_check = pre_treatment_variance>2*pre_treatment_mean;

for j = 1:length(filename.(StrName).BGminSubMean)
    if strains_excluded_check(j,1) == 0
        BGminSubMean_included(j,:) = filename.(StrName).BGminSubMean(j,:);
        Name_included(j,:) = filename.(StrName).Name(j,:);
        yOrf_included(j,:) = filename.(StrName).yOrf(j,:);

    end
end

% Exclude strains that have a variance > 2-fold (pre-treatment)

rows_with_all_zeros = find(all(BGminSubMean_included == 0,2));
BGminSubMean_included(rows_with_all_zeros, : ) = [];

Name_included = Name_included(~cellfun('isempty', Name_included'));
yOrf_included = yOrf_included(~cellfun('isempty', yOrf_included'));

%% Create new struct object

MMS_high_info_included.BGminSubMean_included = BGminSubMean_included;
MMS_high_info_included.Name_included = Name_included;
MMS_high_info_included.yOrf_included = yOrf_included;

%% Calculate median abundance values before treatment and FC

% Calculate median absorbance pre-treatment

pre_treatment_ab = median(MMS_high_info_included.BGminSubMean_included(:,pre_treatment_time_ind),2);

% Calculate fold-change for each timepoint

FC = (MMS_high_info_included.BGminSubMean_included(:,treatment_time_ind))./pre_treatment_ab;

mean_FC = mean(FC,2);

% Add the above into the new struct object

MMS_high_info_included.FC = FC;
MMS_high_info_included.mean_FC = mean_FC;

%% Check maximum absorbance before-treatment and maximum FC

n = 20; % Choose number of genes to rank

maximum_mean_FC = maxk(mean_FC,n);
maximum_median_ab = maxk(pre_treatment_ab,n);

%% Find strains that have highest absorbance before-treatment

[sharedvals_Names_ab,idx_Names_ab] = intersect(pre_treatment_ab, maximum_median_ab,'sorted'); % Names (Ab)
[sharedvals_yOrf_ab,idx_yOrf_ab] = intersect(pre_treatment_ab, maximum_median_ab,'sorted'); % yOrf (Ab)

Names_top20_ab = flip(MMS_high_info_included.Name_included(idx_Names_ab)); % Names (Ab)
yOrf_top20_ab = flip(MMS_high_info_included.yOrf_included(idx_yOrf_ab)); % yOrf (Ab)

%% Find strains that have highest FC

[sharedvals_Names,idx_Names_FC] = intersect(mean_FC, maximum_mean_FC,'sorted'); % Names (FC)
[sharedvals_yOrf,idx_yOrf_FC] = intersect(mean_FC, maximum_mean_FC,'sorted'); % yOrf (FC)

Names_top20_FC = flip(MMS_high_info_included.Name_included(idx_Names_FC)); % Names (FC)
yOrf_top20_FC = flip(MMS_high_info_included.yOrf_included(idx_yOrf_FC)); % yOrf (FC)



%% Find genes with mean FC < 1

housekeeping_gene_finder = find (MMS_high_info_included.mean_FC(idx_Names_ab)<1)

%% Indicate genes' names

ref_genes =(MMS_high_info_included.Name_included(idx_Names_ab))
housekeeping_ref = ref_genes (housekeeping_gene_finder)

%% find pre_treatment abundance

for i = 1:length(housekeeping_ref)
    x = find (strcmp (filename.MMS_high_info.Name,housekeeping_ref(i)))
    pre_treatment_ab_housekeeping(i,1) = pre_treatment_mean(x)
end

for i = 1:length(Names_top20_FC)
    y = find (strcmp (filename.MMS_high_info.Name,Names_top20_FC(i)))
    pre_treatment_ab_genesFC(i,1) = pre_treatment_mean(y)
end