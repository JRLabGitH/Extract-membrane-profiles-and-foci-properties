% Cortex intensity profile extraction from midplane images
% Input: a 60-pixel-wide stripe encompassing the cell membrane the zygote (starting and finishing at the anterior pole).

%% Set variables

% Indicate the channels to analyse
Ch1='false';
Ch2='true'; % will be C3 which is P3
Ch3='true'; % will be C4 which is P2

% Declare and empty
PeakPosAll = struct([]);
DL = struct([]);
DR = struct([]);
D = struct([]);
DC1 = struct([]);
DC2 = struct([]);
DC3 = struct([]);

bdr = struct([]);
cinf1 = struct([]);
cinf2 = struct([]);
cinf3 = struct([]);

ch1M = []; ch2M = []; ch3M =[];
%C1ssh = []; C2ssh = []; C3ssh = [];

% Fitting
%Set channel for inflexion point (InfxC1 InfxC2 or InfxC3). Use in line 260
Infx='InfxC2';

Ch1='false';
Ch2='true'; % will be C3 which is P3
Ch3='true'; % will be C4 which is PK or P2

%For profile representation
maxshift = 0;
extent = 100; % at the moment set extent and extend1 to the same value. Hence normal profiles and extended are the same. %80 %150, it is the pixels you consider at either side of the Infx point
%umpx = 6.60/60; %68.52/512;
smoothInterval = 50; % can be set to 1 if you don't want to smooth the graph

%% Read images and extract A-P profiles to store

% select parent folder to begin analyses and make that the current folder - contains all the images (with all channels) for a treatment
Imgpath = uigetdir;
cd (Imgpath);
F = dir ('*.tif');
info(1).Imgpath = Imgpath;

for f = 1:length(F)     

    fid = F(f).name;
    
    I=imread(fid);
    rowcols = size(I); %extract rows and columns
    
    % Orient images with the outside(top), cortex and cytoplasm(bottom)
            h = fspecial('motion', 9);
            ISm = imfilter(I,h,'Replicate');
            if mean(ISm(1,:)) > mean(ISm(end,:))        
                ISm = flipud(ISm);
                I = flipud(I);
            end
            
    % Extract the cortex
    Project = mean(I,2); % mean of rows to identify where the cortex is
            [PeakVal,PeakPos] = max(Project); 
            PeakPosAll(f).filename = fid;
            PeakPosAll(f).position = PeakPos;
            PeakPos = 33; %activate if the previous step is not good at selecting the cortex (ie. too much external background)
            Cortex = I(PeakPos-3:PeakPos+8, :); %10px on close proximity to the cortex peak
            Cortex = sort(Cortex, 1, 'descend');
            meanCtx = mean(Cortex(1:4,:),1); % mean of the 4 strongest pixels 
    
    % Extract the cytoplasm
            meanCyto = mean(I(end-10:end,:),1); % mean of the 10 bottom pixels 
    
    % Cortex normalisation
            norm_ratio = meanCtx./meanCyto;
            norm_ratio(isinf(norm_ratio)|isnan(norm_ratio))=0; % replaces with zero Inf and NaN
    
    % Cut cortex profiles in half (with end tails) and orient anterior to the left and posterior to the right.
    half_Length = round(length(norm_ratio)/2);
    
    norm_ratiol = cat (2 , norm_ratio(end-150:end) , norm_ratio(1:half_Length+150));
    norm_ratior = fliplr(cat (2 , norm_ratio(half_Length-150:end) , norm_ratio(1:150)));
    
    % Store profiles
    % populate a structure array with normalised half cortical profiles
    % Note that the profiles are of different length and structure array D contains first the left and then the right halves)
    DL(f).filename =cat(2,'L_',fid);
    DL(f).profile = norm_ratiol;   

    DR(f).filename =cat(2,'R_',fid);
    DR(f).profile = norm_ratior;

    D = cat(2, DL, DR);
end


%% Generates a structured array for each channel with every filename linked to a profile

a=1;
b=1;
c=1;
for i=1:length(D)    
   if strfind(D(i).filename,'1-1')
       DC1(a).filename = D(i).filename;
       DC1(a).profile = D(i).profile;
       a=a+1;
   elseif strfind(D(i).filename,'3-1')
       DC2(b).filename = D(i).filename;
       DC2(b).profile = D(i).profile;
       b=b+1;
    elseif strfind(D(i).filename,'4-1')   
       DC3(c).filename = D(i).filename;
       DC3(c).profile = D(i).profile;
       c=c+1;
    end
end

%% CENTRE ALL CHANNELS TO THE INFLEXION POINT OF SELECTED CHANNEL (PAR-3)
    
for i=1:length(DC2) %for i=1:length(DC3)

     if strfind(Ch1,'true')
        C1 = DC1(i).profile;
        filename1 = DC1(i).filename;
     end
     if strfind(Ch2,'true')  
        C2 = DC2(i).profile;
        filename2 = DC2(i).filename;
     end
     if strfind(Ch3,'true')  
        C3 = DC3(i).profile;
        filename3 = DC3(i).filename;
     end

% Obtain point of inflection from Ch1

     if strfind(Ch1,'true')
        xi = 1:length(C1);       
        EFit     =   fittype('a - b * erf(m * (x - c))'); % fitting equation 
        EFitOpts =   fitoptions(EFit);
        
        % set fitting points
        EFitOpts.StartPoint= [0.5 0.5 round(xi(end)/2) 0.1];
        EFitOpts.Upper= [5 5 max(xi) 1];
        EFitOpts.Lower= [0 0 min(xi) 0];
        
        [EFitRes, gof] = fit(xi',C1',EFit,EFitOpts);% retrives the fitted curve and the goodness of fit  
         
        % Generates a structure array to check boundary points 
        bdr(i).a = EFitRes.c; % inflection point
        bdr(i).b = DC1(i).filename;
        
        % Refine: btain the inflexion point of selected curve when only considering points +-200 from the boundary
        boundary = EFitRes.c;        
        k = xi > boundary - 200 & xi < boundary + 200;
        CurveK = C1(k);
        xik = xi(k);
        EFit     =   fittype('a - b * erf(m * (x - c))');
        EFitOpts =   fitoptions(EFit);
    
        EFitOpts.StartPoint= [0.5 0.5 boundary 0.1];
        EFitOpts.Upper= [5 5 boundary+150 1];
        EFitOpts.Lower= [0 0 boundary-150 0];
       
        [EFitRes, gof] = fit(xik',CurveK',EFit,EFitOpts);   
        CurveShift = EFitRes.c;
       
        % Store
        cinf1(i).ptinfx = EFitRes.c; % the refine inflection point
        cinf1(i).filename = DC1(i).filename;
     end

% Obtain point of inflection from Ch2

     if strfind(Ch2,'true')
        xi = 1:length(C2);       
        EFit     =   fittype('a - b * erf(m * (x - c))'); 
        EFitOpts =   fitoptions(EFit);
        
        EFitOpts.StartPoint= [0.5 0.5 round(xi(end)/2) 0.1];
        EFitOpts.Upper= [5 5 max(xi) 1];
        EFitOpts.Lower= [0 0 min(xi) 0];
        
        [EFitRes, gof] = fit(xi',C2',EFit,EFitOpts);% retrives the fitted curve and the goodness of fit
        
        
        % Generates a structure array to check boundary points 
        bdr(i).a = EFitRes.c; %inflection point
        bdr(i).b = DC2(i).filename;
        
        % Refine inflection point
        boundary = EFitRes.c;        
        k = xi > boundary - 200 & xi < boundary + 200;
        CurveK = C2(k);
        xik = xi(k);
        EFit     =   fittype('a - b * erf(m * (x - c))'); 
        EFitOpts =   fitoptions(EFit);
         
        EFitOpts.StartPoint= [0.5 0.5 boundary 0.1];
        EFitOpts.Upper= [5 5 boundary+150 1];
        EFitOpts.Lower= [0 0 boundary-150 0];
        
        [EFitRes, gof] = fit(xik',CurveK',EFit,EFitOpts); 
        CurveShift = EFitRes.c;
         
        % Store
        cinf2(i).ptinfx = EFitRes.c; % the refine inflection point (x value)
        cinf2(i).filename = DC2(i).filename;
     end

% Obtain point of inflection from Ch3

    if strfind(Ch3,'true')
        xi = 1:length(C3);       
        EFit     =   fittype('a - b * erf(m * (x - c))'); 
        EFitOpts =   fitoptions(EFit);
        
        EFitOpts.StartPoint= [0.5 0.5 round(xi(end)/2) 0.1];
        EFitOpts.Upper= [5 5 max(xi) 1];
        EFitOpts.Lower= [0 0 min(xi) 0];
         %ADAPTED FOR POSTERIOR PARS
          %EFitOpts.Upper= [5 5 max(xi) 0];
          %EFitOpts.Lower= [0 0 min(xi) -1];
        
        [EFitRes, gof] = fit(xi',C3',EFit,EFitOpts);% retrives the fitted curve and the goodness of fit
        
        % Generates a structure array to check boundary points 
        bdr(i).a = EFitRes.c;
        bdr(i).b = DC3(i).filename;
        
        % Refine inflection point
        boundary = EFitRes.c;        
        k = xi > boundary - 200 & xi < boundary + 200;
        CurveK = C3(k);
        xik = xi(k);
        EFit     =   fittype('a - b * erf(m * (x - c))'); 
        EFitOpts =   fitoptions(EFit);
         
        EFitOpts.StartPoint= [0.5 0.5 boundary 0.1];
        EFitOpts.Upper= [5 5 boundary+150 1];
        EFitOpts.Lower= [0 0 boundary-150 0];
         %ADAPTED FOR POSTERIOR PARS
          %EFitOpts.Upper= [5 5 boundary+150 0];
          %EFitOpts.Lower= [0 0 boundary-150 -1];
         
        [EFitRes, gof] = fit(xik',CurveK',EFit,EFitOpts);   
        CurveShift = EFitRes.c;
         
        % Store
        cinf3(i).ptinfx = EFitRes.c; % the refine inflection point
        cinf3(i).filename = DC3(i).filename;
    end

 % Sets the Infx point (shift) as 0 in x (centered it in the x axis).
 
     if strcmp(Infx, 'InfxC1') 
         xi = xi-(cinf1(i).ptinfx);
         shift = find(xi>0,1); % rounds up shift
    
     elseif strcmp(Infx, 'InfxC2') 
         xi = xi-(cinf2(i).ptinfx);
         shift = find(xi>0,1);
    
     elseif strcmp(Infx, 'InfxC3') 
         xi = xi-(cinf3(i).ptinfx);
         shift = find(xi>0,1); 
    
     end

    if shift > maxshift
        maxshift = shift;
    end


%For each channel collects profiles in matrixes. Profiles ranging [-extend:extend] with the 0 centered in selected channel inflection
     if strfind(Ch1,'true')
       sC1= smooth(C1,smoothInterval);
       ch1M(i,:) = sC1(shift-extent:shift+extent);% matrix for ploting  - cut profile centred around seleted channel point of inflection and smoothened
       %C1ssh(i,:) = C1(shift-extent:shift+extent); % cut profile and centred around the selected channel inflection point
     end

     if strfind(Ch2,'true')
       sC2 = smooth(C2,smoothInterval);
       ch2M(i,:) = sC2(shift-extent:shift+extent);
       %C2ssh(i,:) = C2(shift-extent:shift+extent); 
     end

     if strfind(Ch3,'true')
        sC3 = smooth(C3,smoothInterval);
        ch3M(i,:) = sC3(shift-extent:shift+extent);
        %C3ssh(i,:) = C3(shift-extent:shift+extent);
     end


end

%% Average profiles from the same embryo 

rowcolsM = size(ch3M);
c = rowcolsM(1)/2;

for i=1:c
    if strfind(Ch1,'true')
        a = ch1M(i,:);
        b = ch1M(i+c,:);
        meanEmb = (a+b)/2; % mean of profiles from same embryos
        ch1Me(i,:) = meanEmb;
     end

     if strfind(Ch2,'true')
        a = ch2M(i,:);
        b = ch2M(i+c,:);
        meanEmb = (a+b)/2; 
        ch2Me(i,:) = meanEmb; 
     end

     if strfind(Ch3,'true')
        a = ch3M(i,:);
        b = ch3M(i+c,:);
        meanEmb = (a+b)/2; 
        ch3Me(i,:) = meanEmb; 
     end
end
%% Adds new data to DC Data structure
for i=1:length(DC2)
    
     if strfind(Ch1,'true')
        DC1(i).ptinfx = cinf1(i).ptinfx; % refined pntinf
        DC1(i).profileFix = ch1M(i,:); % contains profiles cut around point of inflexion and smoothen for plotting
       % DC1(i).profileFixSt = C1ssh(i,:); % contains profiles cut around point of inflexion of selected channel
     end
     if strfind(Ch2,'true') 
        DC2(i).ptinfx = cinf2(i).ptinfx; 
        DC2(i).profileFix = ch2M(i,:); 
       % DC2(i).profileFixSt = C2ssh(i,:); 
     end
     if strfind(Ch3,'true') 
        DC3(i).ptinfx = cinf3(i).ptinfx; 
        DC3(i).profileFix = ch3M(i,:); 
       % DC3(i).profileFixSt = C3ssh(i,:);
     end
end

% Keep only the useful variables
clearvars -except ch1Me ch2Me ch3Me DC1 DC2 DC3