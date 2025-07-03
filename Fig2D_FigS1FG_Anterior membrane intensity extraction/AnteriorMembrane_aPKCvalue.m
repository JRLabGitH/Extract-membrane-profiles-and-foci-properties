clear all

F = dir ('*.tif');
for f = 1:length(F)     
    %Store the image in a matrix where each column is 60px hight (width of
    %selection)
    fid = F(f).name;
    I=imread(fid);

    
    half_Length = round(length(I)/2);

    %Makingnorm ratio L
    firstcut = I(1:60,1:half_Length); 
    secondcut = I(1:60,half_Length:end);
    norm_ratiol = firstcut;
    norm_ratior = fliplr(secondcut)



    
        DL(f).filename =cat(2,'L_',fid);
        DL(f).profile = norm_ratiol;   
  
        DR(f).filename =cat(2,'R_',fid);
        DR(f).profile = norm_ratior;

         D = cat(2, DL, DR);
end


norm_profs = cell(length(D), 1);
anterior_means = nan(length(D),1); % preallocate for means
embryo_names = strings(length(D),1); % preallocate for names

for f = 1:length(D)
    prof = D(f).profile;
    profctx = prof(20:40, :); 
    profctx = sort(profctx, 'descend');
    profctx = mean(profctx(1:10, :));
    profcyto = mean(mean(prof(50:60, end-100:end-50)));
    norm_prof = profctx / profcyto;

    norm_profs{f} = norm_prof;

    % Store embryo (image) name
    embryo_names(f) = D(f).filename;

    % Calculate mean over positions 20–100 (or as much as possible)
    if length(norm_prof) >= 100
        anterior_means(f) = mean(norm_prof(20:100));
    elseif length(norm_prof) >= 20
        anterior_means(f) = mean(norm_prof(20:end));
    else
        anterior_means(f) = NaN; % not enough data
    end
end

% Create a table and save it
EmbryoSummary = table(embryo_names, anterior_means, 'VariableNames', {'EmbryoName', 'AnteriorMean'});
save('EmbryoAnteriorMeans.mat','EmbryoSummary'); % To save as a MATLAB file

