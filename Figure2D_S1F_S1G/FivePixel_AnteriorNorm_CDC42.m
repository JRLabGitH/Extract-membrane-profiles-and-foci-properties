clear all

% Find all .tif files in the current directory
F = dir('*.tif');

for f = 1:length(F)
    % Read each image file
    fid = F(f).name;
    I = imread(fid);
    
    % Determine halfway point of the image
    half_Length = round(size(I, 2) / 2);
    
    % Split the image into two halves and flip the second half
    firstcut = I(1:60, 1:half_Length); 
    secondcut = I(1:60, half_Length:end);
    norm_ratiol = firstcut;
    norm_ratior = fliplr(secondcut);
    
    % Store left and right halves in structures DL and DR
    DL(f).filename = ['L_', fid];
    DL(f).profile = norm_ratiol;
    
    DR(f).filename = ['R_', fid];
    DR(f).profile = norm_ratior;
end

% Concatenate DL and DR into a single structure D
D = [DL, DR];

% Pre-allocate cell arrays and vectors
norm_profs = cell(length(D), 1);
anterior_means = nan(length(D), 1);
embryo_names = strings(length(D), 1);

for f = 1:length(D)
    % Process each profile
    prof = D(f).profile;
    profctx = prof(20:40, :); 
    profctx = sort(profctx, 'descend');
 
    
    profctx = mean(profctx(1:5, :))
    profcyto = mean(mean(prof(57:60, 20:80)))
    

    norm_prof = profctx / profcyto;

    norm_profs{f} = norm_prof;

    % Store embryo (image) name
    embryo_names(f) = D(f).filename;

    % Calculate mean anterior value
    if length(norm_prof) >= 100
        anterior_means(f) = mean(norm_prof(20:100));
    elseif length(norm_prof) >= 20
        anterior_means(f) = mean(norm_prof(20:end));
    else
        anterior_means(f) = NaN; % Not enough data
    end
end

% Create a table and save it
EmbryoSummary = table(embryo_names, anterior_means, ...
    'VariableNames', {'EmbryoName', 'AnteriorMean'});
save('EmbryoAnteriorMeans.mat', 'EmbryoSummary');