strainDir=dir;
allData=struct;
for s=3:length(strainDir)
    cd(strainDir(s).name)
    allData(s).name=strainDir(s).name;
    tiffDir=dir('*.tif');
    allFociProps=[];
    for t=1:length(tiffDir)
        imageData=imEx1(tiffDir(t).name);
        cd('masks')
        maskName=dir(strcat('*',tiffDir(t).name(1:end-4),'*'));
        mask=imEx1(maskName(end).name);
        cd ..
        %     figure;
        %     imshow(imageData(:,:,1),[])
        %     hold on
        [row,col]=find(bwperim(mask));
        frameFociProps=[];
        for tp=1:5
        frame2use=imageData(:,:,tp*10);
        th=multithresh(mat2gray(frame2use),2)
        testMask=imbinarize(mat2gray(frame2use),th(2));
        %figure; imshow(bwareaopen(testMask,3),[])
        bw21=bwareaopen(testMask,3);
        bw3 = imopen(bw21, strel('disk',2) );
        bw4=bw3.*logical(mask);
        % add erosion step
        cc = bwconncomp(bw4,4);
        L_cells = labelmatrix(cc);
        % Add code which skips through frames in 10s e.g. 10,20,30... adds
        % all foci into big matrix, calculates mean? And mean no foci removing when none found?
        figure(1)
        subplot(1,2,1)
        imshow(imageData(:,:,1),[])
        subplot(1,2,2)
        imshow(label2rgb(L_cells,'jet','k','shuffle'));
        hold on
        scatter(col,row)
       % saveas(1,strcat('E:\Data\C elegans\NMY2 spinning disk\Aaron reorganised data\Fig5_Establishment\plots2\',tiffDir(t).name))
        close(1)
        RP=regionprops(L_cells,'Area','circularity');
        fociProps=[];
        for r=1:length(RP)
            fociProps(1,r)=RP(r).Area;
            fociProps(2,r)=RP(r).Circularity;
            fociProps(3,r)=t;

            thisMask=L_cells==r;
            fociPix=uint16(thisMask).*frame2use;
            fociProps(4,r)=mean(fociPix(fociPix>0));
            allMask=L_cells>0;
            backgroundPix=(uint16(logical(mask)-allMask).*frame2use);
            fociProps(5,r)=mean(backgroundPix(backgroundPix>0));
            fociProps(6,r)=fociProps(4,r)/fociProps(5,r);
        end
        
           frameFociProps=cat(2,frameFociProps, fociProps);

        end
        allFociProps=cat(2,allFociProps, fociProps);
        try
            figure;
            KDFplot(frameFociProps(1,:))
            pause
            tiffDir(t).meanFociSize=mean(frameFociProps(1,frameFociProps(1,:)>50));
            tiffDir(t).stdFociSize=std(frameFociProps(1,:));
            tiffDir(t).meanFociIntRatio=mean(frameFociProps(6,:));
            tiffDir(t).stdFociIntRatio=std(frameFociProps(6,:));
            tiffDir(t).Nfoci=length(frameFociProps);
            tiffDir(t).NfociPerFrame=length(frameFociProps)/tp;

        catch
        end
    end
    allData(s).allFociProps=allFociProps;
    T=struct2table(tiffDir);
    filename=strcat('E:\Data\C elegans\NMY2 spinning disk\Aaron reorganised data\Fig5_Establishment\excel files4\',strainDir(s).name,'.xlsx')
    writetable(T,filename)
    allData(s).allFociProps=allFociProps;
    cd ..
end

figure;
barErrorSigJitter({allData(4).allFociProps(1,:),allData(3).allFociProps(1,:),...
    allData(5).allFociProps(1,:),allData(6).allFociProps(1,:),...
    allData(7).allFociProps(1,:),allData(8).allFociProps(1,:)},'t')

xticklabels({allData(4).name,allData(3).name,...
    allData(5).name,allData(6).name,...
    allData(7).name,allData(8).name})


figure;
barErrorSigJitter({allData(4).allFociProps(1,allData(4).allFociProps(1,:)<3000 & allData(4).allFociProps(2,:)>0.5),...
    allData(3).allFociProps(1,allData(3).allFociProps(1,:)<3000 & allData(3).allFociProps(2,:)>0.5),...
    allData(5).allFociProps(1,allData(5).allFociProps(1,:)<3000 & allData(5).allFociProps(2,:)>0.5),...
    allData(6).allFociProps(1,allData(6).allFociProps(1,:)<3000 & allData(6).allFociProps(2,:)>0.5),...
    allData(7).allFociProps(1,allData(7).allFociProps(1,:)<3000 & allData(7).allFociProps(2,:)>0.5),...
    allData(8).allFociProps(1,allData(8).allFociProps(1,:)<3000 & allData(8).allFociProps(2,:)>0.5)},'t')

xticklabels({allData(4).name,allData(3).name,...
    allData(5).name,allData(6).name,...
    allData(7).name,allData(8).name})

figure;
barErrorSigJitter({allData(4).allFociProps(6,:),allData(3).allFociProps(6,:),...
    allData(5).allFociProps(6,:),allData(6).allFociProps(6,:),...
    allData(7).allFociProps(6,:),allData(8).allFociProps(6,:)},'t')

xticklabels({allData(4).name,allData(3).name,...
    allData(5).name,allData(6).name,...
    allData(7).name,allData(8).name})