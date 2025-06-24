%% In the command windown you can run this function by running "[allMeans,lowerBound, upperBound, Ydata, Table] =bootstrapProfileFits(profileDataNAME, 10, 1.8)"
% profileData could be ch1Me, ch2Me, ch3Me. RefY is the Y value where you want to measure profile distance

%function [allMeans, lowerBound, upperBound, Ydata, Table] = PlotProfBootstrap(profileData, niterations, RefY)
function [PrsmCurveK, PrsmSD, PrsmxikUm, PrsmNembM, Table] = PlotProfBootstrap(profileData, niterations, RefY)

%% from a matrix of profiles (a profile per row) calculate the mean and SEM 
umpx = 6.60/60;
L= size(profileData,2);
Half = round((L/2)-1);

Nembryos=size(profileData,1);
xik=1:size(profileData,2);
CurveK=mean(profileData,1);

for p=1:size(profileData,2)
    lowerBoundSEM(p)= mean(profileData(:,p)) - std(profileData(:,p))/Nembryos^0.5;
    upperBoundSEM(p)= mean(profileData(:,p)) + std(profileData(:,p))/Nembryos^0.5;
    SD(p)= std(profileData(:,p));
end

% For prism plotting
PrsmCurveK = CurveK';
PrsmSD = SD';
Prsmxik = xik - Half;
PrsmxikUm = Prsmxik * umpx;
PrsmxikUm = PrsmxikUm';
NembM = xik;
NembM(NembM>0) = Nembryos;
PrsmNembM = NembM';

%% Fit curve with the erf function to the data
EFit     =   fittype('a - b * erf(m * (x - c))'); % fitting equation 
EFitOpts =   fitoptions(EFit);
EFitOpts.StartPoint= [0.5 0.5 100 0.1];
EFitOpts.Upper= [Inf Inf Inf Inf];
EFitOpts.Lower= [-Inf -Inf -Inf -Inf];
[EFitRes, gof] = fit(xik',CurveK',EFit,EFitOpts); 
%plot(EFitRes, 'predobs'); % predobs contains the prediction bounds
Ydata = EFitRes(xik);

% Extract the X value for the Y position where the distance between the
% profiles wants to be calculated (RefY)
ftogetx = @(x) RefY - EFitRes(x); % for every x value calculate the y of profile and subtract it by RefY
xfound = fzero(ftogetx, 0); % find the x values where ftogetx is 0, this matches the X in profile at the RefY


%% Fit curve with the erf function to obtain the shift needed to centre the point of inflection
k = xik >40 & xik <160; %to exclude 40 points at either end (x >40 & x <160)
xik = xik(k);
CurveK = CurveK(k);

EFit     =   fittype('a - b * erf(m * (x - c))'); % fitting equation 
EFitOpts =   fitoptions(EFit);
EFitOpts.StartPoint= [0.5 0.5 100 0.1];
EFitOpts.Upper= [Inf Inf Inf Inf];
EFitOpts.Lower= [-Inf -Inf -Inf -Inf];
[EFitRes, gof] = fit(xik',CurveK',EFit,EFitOpts); 
%plot(EFitRes, 'predobs'); % predobs contains the prediction bounds
pntInf = EFitRes.c; % to determine if we need to readjust the the PAR-3 mean curve to 0 (should be a value of a 100) - fitting should be done when reomving 40points at either end -
shift = Half-pntInf;
shift = round(shift);

%reset for later
xik=1:size(profileData,2);
CurveK=mean(profileData,1);

%% Bootstrapping
% This section generates new profile matrices (newProfile) by randomly sampling profiles from the original data with replacement until the matrix reaches the same size as the original.
% For each resampled matrix, it calculates the mean profile and fits a curve. This process is repeated for the number of iterations specified, and all fitted curves are stored in allMeans.

%figure;
allMeans=zeros(niterations, size(profileData,2));
for n=1:niterations
    randInd=randi([1, Nembryos],1,Nembryos);
    newProfile=profileData;
    for p=1:Nembryos
        newProfile(p,:)=profileData(randInd(p),:);
    end
    CurveK=mean(newProfile,1);
    [EFitResBoot] = fit(xik',CurveK',EFit,EFitOpts);
    %allMeansNF(n,:)=CurveK;
    allMeans(n,:)=EFitResBoot(xik);
    %plot(allMeans(n,:))
    %hold on
end

for p=1:size(profileData,2)
    lowerBound(p)=prctile(allMeans(:,p),5);
    upperBound(p)=prctile(allMeans(:,p),95);
end

% Extract the X error (95%CI from boostrapfit data) 
ErrorCurveL=lowerBound';
[~,minInd]=min(abs(ErrorCurveL(:,1)-EFitRes(xfound)));
CIfoundL=xfound-xik(minInd);

ErrorCurveU=upperBound';
[~,minInd]=min(abs(ErrorCurveU(:,1)-EFitRes(xfound)));
CIfoundU=xfound-xik(minInd);
CIfoundU=abs(CIfoundU);

CIfound = (CIfoundL+CIfoundU)/2; % 95%CI used to add in quadrature to calculate the error for the distance between to profiles at the given RefY

%% Generate output

varNames = {'Y','xfound','95CIBootStrp','pntInflec', 'shift'};
Table = table (RefY, xfound, CIfound, pntInf, shift, 'VariableNames',varNames);

figure; 
errorbar(xik,mean(profileData,1),std(profileData,0,1)/Nembryos^0.5)
hold on
plot(EFitRes, 'predobs')
hold on
plot(mean(allMeans));
hold on
plot(1:size(profileData,2),upperBound);
plot(1:size(profileData,2),lowerBound);

legend('Mean data','fitted curve to mean data','95% CI on mean fit upper','95% CI on mean fit lower', 'bootstrap mean', '95% CI bootstrap upper', '95% CI bootstrap lower')

end
