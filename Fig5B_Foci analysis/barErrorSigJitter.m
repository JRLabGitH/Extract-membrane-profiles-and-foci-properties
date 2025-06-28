function p=barErrorSigJitter(data,test2run)
%data is a cell array of different raw data for means to be calculated from
%test2run=t for 2 sample t test or w for Wilcoxon rank sum test (good for
%non-normal distributions)
for t=1:size(data,2)
    means(t)=nanmean(data{t});
    SD(t)=nanstd(data{t});
    N(t)=length(data{t});
    if t==1
        ttest2val(t)=0;
        ttest2txt{t}=' ';
    else
        if test2run=='t';
            [ttest2val(t),p(t)]=ttest2(data{1},data{t});
        elseif test2run=='w'
            
            [p(t),ttest2val(t)] = ranksum(data{1},data{t});
        elseif test2run=='b'
            p(t) = brunner_munzel(data{1}, data{t});
            ttest2val(t)=p(t)<0.05;
        end
        if ttest2val(t)==1
            ttest2txt{t}='*';
        else
            ttest2txt{t}=' ';
        end
    end
end


bar(means,'FaceColor','w')

%boxplot(means,'Colors','k')
hold on
% for t=1:size(data,2)
%     jitterVar=0.4*rand(size(data{t}))-0.2+t;
% scatter(jitterVar,data{t},0.05,'r.')
% end

for t=1:size(data,2)
    jitterVar=0.4*rand(size(data{t}))-0.2+t;
    dataToUse=data{t};
scatter(jitterVar(1:end),dataToUse(1:end),1,'k.')
end

bar(means,'FaceAlpha',0.5,'FaceColor',[.7 .7 .7])
 errorbar(1:length(means),means,SD./N.^0.5,'k','LineStyle','none')
%errorbar(1:length(means),means,SD,'k','LineStyle','none')

%errorbar(1:length(means),means,SD./N.^0.5,'k')
text((1:length(means))-0.025,(means+SD./N.^0.5)*1.2,ttest2txt,'FontSize',20)

end