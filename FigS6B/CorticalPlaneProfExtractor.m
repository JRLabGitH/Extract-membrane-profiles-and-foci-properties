%% Assemble data structure (cortex) with matrices from all images

%Find cortical bands to analyse
loc= "/Users/josana/Desktop/FigS5B_CTXanalyses_upld/Images/aPKC cortical bands";
tr_dir = dir(loc);
tr_dir(strncmp({tr_dir.name}, '.', 1)) = []; % skip files and dir starting with '.'
sc = 1;

%Organizes image matrices into a data structure based on the folder hierarchy
%loop through each condition folder
for condition = 1:length(tr_dir)
    
    %new path into condition folder
    ori_path = strcat(loc, '/', tr_dir(condition).name);
    nu_dir = dir(ori_path);
    nu_dir(strncmp({nu_dir.name}, '.', 1)) = []; % skip files and dir starting with '.'
    
    %loop through genotypes
    %for gen = 3:length(nu_dir)
    
        %find and loop through image folders
        %fin_path = strcat(ori_path, '/', nu_dir(gen).name);
        %fin_dir = dir(fin_path);
        
        %loop through images
        for ima = 1:length(nu_dir)
            %filepath
            im_path = strcat(ori_path,'/', nu_dir(ima).name);
            %im_path = strcat(fin_path,'/', fin_dir(ima).name);
            
            %read image into structure for analysis
            IM = imread(im_path);
           
            %read into structure
            cortex(sc).Strain = tr_dir(condition).name; %strain/genotype
            %cortex(sc).RNAi = nu_dir(gen).name; %RNAi
            cortex(sc).Embryo = nu_dir(ima).name; %Image/embryo
            cortex(sc).Image = IM;
            
            sc = sc+1;
        end
    %end
end

%% Extract AP Intensity Profile from cortex band

ROI_width = 50; %size of averaging ROI

for n =1:length(cortex)
    
    IM = cortex(n).Image;
    
    cl = length(IM); %cortex length
    
    bc_int=[]; %initialise array
    
    % Move a window of defined width (ROI_width) across the image band one pixel at a time.
    % At each position, calculate the average intensity of the pixels within the window.
    
    for n1 = 1:cl-ROI_width
        ROI = IM(:,n1:n1+ROI_width-1); %define ROI
        bc_int(n1) = mean(mean(ROI, 1)); % obtain meain intensity for ROI
    end  

    %Store AP intensity profile
    cortex(n).AP_intensity = bc_int;
    
    %plot
%         figure
%         plot(bc_int)
    
  
    %bin AP intensity profile in percentiles - (allows co-plotting  domain
    %length as percentage length of embryos) 
    pc = length(bc_int)/100; %one percent
    percent = [];
    
    for n1 = 1:100
        start_i = (n1-1)*pc +1;
        end_i = n1*pc;
        percent(n1,1) = start_i; %start index of bin
        percent(n1,2) = end_i;   %end index of bin
        %average AP profile over the bin
        percent(n1,3) = mean(bc_int(floor(start_i):floor(end_i)));
    end
    
    %store percentiled cortex in cortex
    cortex(n).percentile_cortex = percent;
    
    %plot percentiled cortex profile
%     figure
%     plot(percent(:,3))
%     hold on
end


%% cortex data extractor into s2

Strains = {'CDC42WT', 'CDC42S71A', 'CDC42S71E'};
RNAis = {'utr'};%, 'utr 50', 'chin1 50'};
channels = { 'C2'}; %'C1',

% Group the percentiled cortex per experimental condition
hc = 0; s2 = {};
for str = 1:length(Strains)
    strn = Strains{str};
    
    for kd = 1:length(RNAis)
        rn = RNAis{kd};
        
        for ch = 1:length(channels)
            cha = channels{ch};
            
            hc = hc+1;
            c=0; s1 = {};
            for n = 1:length(cortex)
                
                strain = cortex(n).Strain;
                %RNAi = cortex(n).Strain;
                file = cortex(n).Embryo;
                
                if strcmp(strain,strn) && contains(file,cha)
                %if strcmp(strain,strn) &&  strcmp(RNAi,rn) && contains(file,cha)
                    
                    c=c+1;
                    s1{c} = cortex(n).percentile_cortex;
                end
                
            end
            
            %string = strcat(Strains{str},'_', RNAis{kd}, '_', channels{ch});
            string = strcat(Strains{str},'_', channels{ch});
            s2(hc).id = string;
            s2(hc).percent = s1;
        end
    end
end

%build s2
for n = 1:length(s2) % for each category
    p_int = []; %initialise array
    
    for n1 = 1:length(s2(n).percent) %for each profile
        p_int(:,n1) = s2(n).percent{n1}(:,3);
    end
    
    av_int = mean(p_int,2); %calculate AP prof mean for each condition
    av_int = av_int'; %flip to row
    
    for n1 = 1:100 % for each percent
        sd_int(n1) = std(p_int(n1,:));
    end
    
    %store in s2

    s2(n).int_percent_table = p_int;
    s2(n).int_mean_percent = av_int;
    s2(n).int_std_dev = sd_int;
end

%Standard error of the Mean (SEM)
for n = 1:length(s2)
    
    n_embryos = length(s2(n).percent);
   
    std_dev_int = s2(n).int_std_dev;
    int_SEM = std_dev_int/sqrt(n_embryos);
    s2(n).int_SEM = int_SEM;
end

%% Set HiLo ratio for plotting AP Intensity profiles

for n = 1:length(s2)
    
    int = s2(n).int_percent_table;
    
    embryo_n = length(int(1,:));
    
    hilo_int =[];
    
    for n1 = 1:embryo_n
        
        
        int_data = int(:,n1); %collect embryo data
        
        %apply ratio to int data
        int_sort = sort(int_data);
        int_av_lo = mean(int_sort(1:10)); int_av_hi = mean(int_sort(91:100));
        ratio_int = int_data - int_av_lo;
        norm_int = ratio_int/(int_av_hi-int_av_lo);
        
        hilo_int(:,n1) = norm_int;
        av_hilo_int = mean(hilo_int, 2);
        std_int = std(hilo_int,1,2);
        sem_int = std_int/sqrt(embryo_n);
    end 
    
    s2(n).hilo_int = hilo_int;
    s2(n).mean_hilo_int = av_hilo_int;
    s2(n).hilo_int_std = std_int;
    s2(n).hilo_int_sem = sem_int;
end


%% PLOTTING

wt = [0 0.4470 0.7410];
S71A = 	[1 0 1];
S71E = [0.9290 0.6940 0.1250];

%SETTINGS
% 1 for standard dev  2 for SEM
error = 2;

%collect data


s1data = s2(1).mean_hilo_int;
s2data = s2(2).mean_hilo_int;
s3data = s2(3).mean_hilo_int;

if error == 1
    er1 = s2(1).hilo_int_std;
    er2 = s2(2).hilo_int_std;
    er3 = s2(3).hilo_int_std;
    
else
    er1 = s2(1).hilo_int_sem;
    er2 = s2(2).hilo_int_sem;
    er3 = s2(3).hilo_int_sem;
end




figure
%wt
plot(s1data, 'Color', wt, 'LineWidth', 2)
hold on
plot(s1data+er1,'Color', wt,'LineWidth', 0.5)
hold on
plot(s1data-er1,'Color', wt, 'LineWidth', 0.5)
hold on
%S71A
plot(s2data, 'Color', S71A, 'LineWidth', 2)
hold on
plot(s2data+er2,'Color', S71A,'LineWidth', 0.5)
hold on
plot(s2data-er2,'Color', S71A, 'LineWidth', 0.5)
hold on
%S71E
plot(s3data, 'Color', S71E, 'LineWidth', 2)
hold on
plot(s3data+er3,'Color', S71E,'LineWidth', 0.5)
hold on
plot(s3data-er3,'Color', S71E, 'LineWidth', 0.5)
hold on

% Keep only the useful variables
clearvars -except cortex s2
