
%% Store image matrices in the structured array "cortex" for downstream analysis
% Define location of data
%loc = "C:\Users\brook\OneDrive\Desktop\par6 50_script_test_can_delete"; %pc
loc = "/Users/josana/Desktop/aliciaJR"; %mac


% Get condition folders (A, E, WT)
tr_dir = dir(loc);
tr_dir = tr_dir([tr_dir.isdir]); % Keep only directories
tr_dir = tr_dir(~ismember({tr_dir.name}, {'.', '..'})); % Remove '.' and '..'

sc = 1; % Counter for struct

% Loop through condition folders (A, E, WT)
for condition = 1:length(tr_dir)
    
    % New path into condition folder (A, E, or WT)
    condition_path = fullfile(loc, tr_dir(condition).name);
    
    % Find all TIF images in this folder
    image_files = dir(fullfile(condition_path, '*.tif'));
    
    % Debugging: Check if images exist
    disp(['Processing folder: ', condition_path]);
    disp(['Number of images found: ', num2str(length(image_files))]);
    
    if isempty(image_files)
        warning(['No images found in: ', condition_path]);
        continue;
    end
    
    % Loop through images
    for img_idx = 1:length(image_files)
        % File path
        img_path = fullfile(condition_path, image_files(img_idx).name);
        
        % Ensure file exists before reading
        if exist(img_path, 'file') ~= 2
            warning(['File not found: ', img_path]);
            continue;
        end
        
        % Read image into structure for analysis
        IM = imread(img_path);
       
        % Store data in struct
        cortex(sc).RNAi = tr_dir(condition).name; % RNAi condition (A, E, WT)
        cortex(sc).Strain = "Default"; % Default genotype (since there's no genotype subfolder)
        cortex(sc).Embryo = image_files(img_idx).name; % Image filename
        cortex(sc).Image = IM; % Store image
        
        sc = sc + 1; % Increment struct index
    end
end

%% Extract Membrane crossing profile and Peak membrane values at anterior and posterior

for n =1:length(cortex)
    
    IM = cortex(n).Image;
    
    cl = length(IM); %cortex length
    
    %anterior ROIs
    ant1 = IM(:,60:120);
    ant2 = IM(:,cl-120:cl-60);

    %posterior ROI
    post = IM(:,(floor(cl/2)-60):(floor(cl/2)+60));
    
    %local cytoplasm
    lb_ant1 = mean(mean(ant1(50:60,:)));
    lb_ant2 = mean(mean(ant2(50:60,:)));
    lb_post = mean(mean(post(50:60,:)));
    
    %background correction
    bc_ant1 = ant1/lb_ant1;
    bc_ant2 = ant2/lb_ant2;
    bc_post = post/lb_post;
    
    %Find profile
    ant_prof1 = mean(bc_ant1,2);
    ant_prof2 = mean(bc_ant2,2);
    %combine anterior profiles
    ant = [ant_prof1, ant_prof2];
    ant_av = mean(ant,2);
    
    pos_prof = mean(bc_post,2);
    
    %find cortex
    %look for first peak around 30 pixels
    peakant1 = islocalmax(ant_av(29:40));
    peakpost = islocalmax(pos_prof(29:40));
    
    ant1_i = 0; post_i = 0;
    
    %check for multiple peaks
    if sum(peakant1) == 1
        
        for n1 = 1:length(peakant1)
            if peakant1(n1) == 1
                ant1_i=28+n1;
                break
            end
        end
        
    else
        %find the peaks and their values
        Ay = ant_av(30:40);
        A = Ay(peakant1 == 1);
        %see if the first peak is in the background
        for n2 = 2:length(A)
            if A(n2) - 0.2 > A(n2-1)
                As = cumsum(peakant1);
                %identify index of n2 peak in As
                for n3 = 1:length(As)
                    if As(n3) == n2
                        ant1_i=28+n3;
                        break
                    end
                end
            end
        end
        
        if ant1_i == 0
            
            for n1 = 1:length(peakant1)
                if peakant1(n1) == 1
                    ant1_i=28+n1;
                    break
                end
            end
        end
        
        if ant1_i == 0
            for n1 = 1:length(Ay)-1
                if Ay(n1+1) < Ay(n1)+0.01
                    ant1_i=28+n1;
                    break
                end
            end
            
        end
    end
        
        
    if sum(peakpost) == 1
        
        for n1 = 1:length(peakpost)
            if peakpost(n1) == 1
                post_i=28+n1;
                break
            end
        end
        
    else
        %find the peaks and their values
        Cy = pos_prof(30:40);
        C = Cy(peakpost == 1);
        %see if the first peak is in the background
        for n2 = 2:length(C)
            if C(n2) - 0.2 > C(n2-1)
                Cs = cumsum(peakpost);
                %identify index of n2 peak in Cs
                for n3 = 1:length(Cs)
                    if Cs(n3) == n2
                        post_i=28+n3;
                        break
                    end
                end
            end
        end
            
        if post_i == 0
            
            for n1 = 1:length(peakpost)
                if peakpost(n1) == 1
                    post_i=28+n1;
                    break
                end
            end
        end
        
        if post_i == 0
            for n1 = 1:length(Cy)-1
                if Cy(n1+1) < Cy(n1)+0.01
                    post_i=28+n1;
                    break
                end
            end
        end
    end

   ppeaki(n,:) = [ant1_i, post_i];
   peakv(n,:) = [ant_av(ant1_i), pos_prof(post_i)];

    % % Plot profile and peak selection
    % figure 
    % plot(ant_av)
    % hold on
    % scatter(ant1_i, ant_av(ant1_i),'k')
    % hold on
    % plot(pos_prof)
    % hold on
    % scatter(post_i, pos_prof(post_i),'k')

    %Save pixel ROIs in cortex
    cortex(n).Anterior1_raw = ant1;
    cortex(n).Anterior1_bkgrnd = bc_ant1;
    cortex(n).Anterior2_raw = ant2;
    cortex(n).Anterior2_bkgrnd = bc_ant2;
    cortex(n).Posterior_raw = post;
    cortex(n).Posterior_bkgrnd = bc_post;
    
    %Save analysis in cortex
   cortex(n).anterior_peak = peakv(n,1);
   cortex(n).anterior_peak_index = ant1_i;
   cortex(n).posterior_peak = peakv(n,2);
   cortex(n).posterior_peak_index = post_i;
   cortex(n).anterior_profile = ant_av;
   cortex(n).posterior_profile = pos_prof;
end


%% Centre profiles by their peaks and average by condition for profile plotting
ant_max = max(ppeaki(:,1));
ant_min = min(ppeaki(:,1));
ant_range = ant_max - ant_min;

pos_max = max(ppeaki(:,2));
pos_min = min(ppeaki(:,2));
pos_range = pos_max - pos_min;
    
ant_prof_centre = nan(length(ppeaki), ant_range+length(ant_av));
pos_prof_centre = nan(length(ppeaki), pos_range+length(pos_prof));

for n = 1:length(ppeaki)
    
   ant_i = ppeaki(n,1);
   gap = (ant_max - ant_i)+1;
   ant_prof_centre(n,gap:gap+length(cortex(n).anterior_profile)-1) = cortex(n).anterior_profile;

    
   pos_i = ppeaki(n,2);
   gap = (pos_max - pos_i)+1;
   pos_prof_centre(n,gap:gap+length(cortex(n).posterior_profile)-1) = cortex(n).posterior_profile;
    
end

for n= 1: length (cortex)
    cortex(n).ant_prof_centre = ant_prof_centre(n, :);
end

for n1=1:length (tr_dir)
    cond=tr_dir(n1).name;
    s=1;
    for n2=1:length (cortex)
        if cond == cortex(n2).RNAi
            prfCent(n1).cond = cond;
            prfCent(n1).matprf(s,:)= cortex(n2).ant_prof_centre;
            s=s+1;
        end
    end
    prfCent(n1).meanPrf = (mean (prfCent(n1).matprf))';
    prfCent(n1).stdPrf = (std (prfCent(n1).matprf))';
    prfCent(n1).countPrf = height (prfCent(n1).matprf);
    prfCent(n1).peakpos = ant_max;

end
         
% Keep only the useful variables
clearvars -except cortex prfCent peakv ant_prof_centre