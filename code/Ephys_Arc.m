%% EphysAnal_Arc

%previously VisCheckKilosort
%After running preproccesing routine writen by N.Karalis 


clear all
addpath(genpath('...\karaniko\Code\data_suite2'))

load('...\data\Ephys\ArcAll.mat') 

%% Z-scored
for ii=1:length(PSTH_opto20Hz2All)
    temp = PSTH_opto20Hz2All{ii};
    mean2sub = mean(mean(temp(:,1:30)));
    std2dev = std(mean(temp(:,1:60)));
    if std2dev == 0
        PSTH_opto20Hz2Allz{ii,1} = temp*0;
    else
        PSTH_opto20Hz2Allz{ii,1}=(temp - mean2sub)/std2dev;
    end
    PSTH_opto20Hz2Allz_ave(ii,:)=mean(PSTH_opto20Hz2Allz{ii,1},1);
end

for ii=1:length(PSTH_CS2All)
    temp = PSTH_CS2All{ii};
    mean2sub = mean(mean(temp(:,1:30)));
    std2dev = std(mean(temp(:,1:60)));
    if std2dev == 0
        PSTH_CS2Allz{ii,1} = temp*0;
    else
        PSTH_CS2Allz{ii,1}=(temp - mean2sub)/std2dev;
    end
    PSTH_CS2Allz_ave(ii,:)=mean(PSTH_CS2Allz{ii,1},1);
end

%% Sorting based on z-score and visualize
%CS
for ii=1:size(PSTH_CS2All,1)
    CS_diff(ii) = mean(PSTH_CS2Allz_ave(ii,31:40)) - mean(PSTH_CS2Allz_ave(ii,21:30));
    CS_sig(ii) = mean(PSTH_CS2Allz_ave(ii,31:40));
end

[AA_CS,CS_sort]=sort(CS_diff,'descend');

for ii=1:size(PSTH_opto20Hz2All,1)
    temp=PSTH_opto20Hz2All{ii};
    Opto_diff(ii) = mean(PSTH_opto20Hz2Allz_ave(ii,31:40)) - mean(PSTH_opto20Hz2Allz_ave(ii,21:30));
    Opto_sig(ii) = mean(PSTH_opto20Hz2Allz_ave(ii,31:40));
end

[AA_Opto,Opto_sort]=sort(Opto_diff,'descend');

figure,
set(gcf,'position',[100,100,400,400])
subplot(5,2,[1 3 5 7]),imagesc(PSTH_CS2Allz_ave(CS_sort,:),[-3 4])
subplot(5,2,[2 4 6 8]),imagesc(PSTH_opto20Hz2Allz_ave(Opto_sort,:),[-3 4])
subplot(5,2,9),

plot(cumsum(histc(CS_sig,[-2.5:0.01:2.5]))/length(CS_sig),'m'),hold on
plot(cumsum(histc(Opto_sig,[-2.5:0.01:2.5]))/length(Opto_sig),'b'),hold on
xticks([0 250 500])
xticklabels({'-2.5','0','2.5'})
xlim([0 500])

[p,m]=signrank(CS_sig, Opto_sig);
title(['sign-rank p = ' num2str(p)])
%subplot(5,2,10)
%plot([1,2],[CS_sig, Opto_sig],'-',Color=[0.5 0.5 0.5 0.1])
%xlim([0.5 2.5])

figure,
set(gcf,'position',[100,100,400,400])
subplot(5,2,[1 3 5 7]),imagesc(PSTH_CS2Allz_ave(CS_sort,:),[-3 4])
subplot(5,2,[2 4 6 8]),imagesc(PSTH_opto20Hz2Allz_ave(CS_sort,:),[-3 4])

%% rank correlation boot strapping

original_sequence = CS_sig; 
n_bootstrap_samples = 100000; 

bootstrap_correlations = zeros(1, n_bootstrap_samples);
for i = 1:n_bootstrap_samples

    shuffled_sequence = original_sequence(randperm(length(original_sequence)));

    correlation = corr(original_sequence', shuffled_sequence', 'type', 'Spearman');
    
    bootstrap_correlations(i) = correlation;
end

observed_correlation = corr(original_sequence', Opto_sig', 'type', 'Spearman');
%p_value = sum(abs(bootstrap_correlations) >= abs(observed_correlation)) / n_bootstrap_samples;
p_value = sum(bootstrap_correlations >= observed_correlation) / n_bootstrap_samples;

figure;
set(gcf,'position',[100,100,600,400])
histogram(bootstrap_correlations, 'Normalization', 'probability','FaceColor', [0.7 0.7 0.7]);hold on
plot([observed_correlation observed_correlation],[0 0.04],'--')
xlabel('Spearman R');
ylabel('observed frequency');

title( ['bootstrapP:' num2str(p_value)]);
%% Sorting based on response amplitudes
%CS
for ii=1:size(PSTH_CS2All,1)
    temp=PSTH_CS2All{ii};
    CS_diff(ii) = mean(mean(temp(:,31:60))) - mean(mean(temp(:,21:30)));
end

[AA,CS_sort]=sort(CS_diff,'descend');

%Opto
for ii=1:size(PSTH_opto20Hz2All,1)
    temp=PSTH_opto20Hz2All{ii};
    Opto_diff(ii) = mean(mean(temp(:,31:60))) - mean(mean(temp(:,21:30)));
end

[AA,Opto_sort]=sort(Opto_diff,'descend');

%Optotag
for ii=1:size(PSTH_optotag2All,1)
    temp=PSTH_optotag2All{ii};
    Optotag_diff(ii) = mean(temp(:,21))/mean(mean(temp)); %10ms bin
end

[AA,Optotag_sort]=sort(Optotag_diff,'descend');

%% All PSTH checking for selection
%{
for ii=1:length(PSTH_optotag2All)
%for ii=1:5
figure
set(gcf,'position',[100,100,800,600])
temp=PSTH_optotag1All{ii};

subplot(2,2,1)
out = RasterRatio(temp,10);
imshow(imcomplement(out)),hold on

plot([600 600], [0 1000], '--', Color=[0.8 0.5 0.5])
title(['Unit:' num2str(ii)])

temp=PSTH_optotag2All{ii};
subplot(2,2,3),
bar(0.5:39.5, mean(temp,1),1)
ylabel('FR(Hz)')

temp=PSTH_opto20Hz1All{ii};

subplot(2,2,2)
out = RasterRatio(temp,200);
imshow(imcomplement(out)),hold on

plot([600 600], [0 1000], '--', Color=[0.8 0.5 0.5])
title(['Unit:' num2str(ii)])

temp=PSTH_opto20Hz2All{ii};
subplot(2,2,4),
bar(0.5:59.5, mean(temp,1),1)
ylabel('FR(Hz)')

end
%}


%% Vis
%{
for kk=1:length(COI)
    ii=COI(kk);
    %ii=CS_sort(kk);
    %ii=Opto_sort(kk);
    %ii=Optotag_sort(kk);

figure,
set(gcf,'position',[100,100,1800,400])

subplot(2,3,1)
%Square_coloring([200 210], [0.5 0.5 1])
out = RasterRatio(PSTH_optotag1All{ii},1);
imshow(imcomplement(PSTH_optotag1All{ii}(:,101:300))),hold on

plot([100 100], [0 400], '--', Color=[0.8 0.8 0.8])
plot([110 110], [0 400], '--', Color=[0.5 0.5 0.8])
title(['CS  Unit:' num2str(ii)])

subplot(2,3,4)
bar(0.5:19.5,mean(PSTH_optotag2All{ii}(:,11:30),1)*(1000/PSTHparam.optotag_bin2),1,'FaceColor',[0.5 0.5 0.5])
ylabel('FR(Hz)')

subplot(2,3,2)
out = RasterRatio(PSTH_CS1All{ii},200);
imshow(imcomplement(out)),hold on

plot([600 600], [0 1000], '--', Color=[0.8 0.5 0.5])
title(['CS  Unit:' num2str(ii)])

subplot(2,3,5),
bar(0.5:59.5, mean(PSTH_CS2All{ii},1),1,'FaceColor',[0.5 0.5 0.5])
ylabel('FR(Hz)')

%subplot(2,4,3)

%out = RasterRatio(PSTH_opto4Hz1All{ii},100);
%imshow(imcomplement(out)),hold on

%plot([600 600], [0 1000], '--', Color=[0.8 0.5 0.5])
%title(['Opto4Hz  Unit:' num2str(ii)])

%subplot(2,4,7)
%bar(0.5:59.5,mean(PSTH_opto4Hz2All{ii},1),1,'FaceColor',[0.5 0.5 0.5])
%ylabel('FR(Hz)')

subplot(2,3,3)
out = RasterRatio(PSTH_opto20Hz1All{ii},100);
imshow(imcomplement(out)),hold on

plot([600 600], [0 1000], '--', Color=[0.8 0.5 0.5])
title(['Opto20Hz  Unit:' num2str(ii)])

subplot(2,3,6)
bar(0.5:59.5,mean(PSTH_opto20Hz2All{ii},1),1,'FaceColor',[0.5 0.5 0.5])
ylabel('FR(Hz)')

end
%}


%% tone / opto responses

for ii = 1:size(PSTH_CS2All)
    CS_base = mean(PSTH_CS2All{ii}(:,1:29));
    CS_sig = mean(PSTH_CS2All{ii}(:,31:40));
    Opto_base = mean(PSTH_opto20Hz2All{ii}(:,1:29));
    Opto_sig = mean(PSTH_opto20Hz2All{ii}(:,31:40));

    CSp(ii) = ranksum(CS_base, CS_sig);
    Optop(ii) = ranksum(Opto_base, Opto_sig);
    CSsign(ii) = (mean(CS_sig) - mean(CS_base))/abs(mean(CS_sig) - mean(CS_base));
    Optosign(ii) = (mean(Opto_sig) - mean(Opto_base))/abs(mean(Opto_sig) - mean(Opto_base));
    CSsignp(ii) = CSp(ii)*CSsign(ii);
    Optosignp(ii) = Optop(ii)*Optosign(ii);
end

%%
Opto_posi = find(Optosignp > 0 & Optosignp < 0.01);
Opto_nega = find(Optosignp < 0 & Optosignp > -0.01);
CS_posi = find(CSsignp > 0 & CSsignp < 0.01);
CS_nega = find(CSsignp < 0 & CSsignp > -0.01);

display(['CS_posi:' num2str(length(CS_posi)/length(PSTH_CS2All))])
display(['CS_nega:' num2str(length(CS_nega)/length(PSTH_CS2All))])

display(['Opto_posi:' num2str(length(Opto_posi)/length(PSTH_CS2All))])
display(['Opto_nega:' num2str(length(Opto_nega)/length(PSTH_CS2All))])


%%
PSTH_CS2All_CS_posi=[];
PSTH_CS2All_CS_nega=[];
PSTH_opto20Hz2All_Opto_posi=[];
PSTH_opto20Hz2All_Opto_nega=[];

figure,
set(gcf,'position',[100,100,800,600])
subplot(2,2,1)
for ii=1:length(CS_posi)
    PSTH_CS2All_CS_posi = [PSTH_CS2All_CS_posi; mean(PSTH_CS2Allz{CS_posi(ii)})]; %z=scored
end

for ii=1:length(CS_nega)
    PSTH_CS2All_CS_nega = [PSTH_CS2All_CS_nega; mean(PSTH_CS2Allz{CS_nega(ii)})];
end

SEM = std(PSTH_CS2All_CS_posi)/sqrt(length(PSTH_CS2All_CS_posi));
%SEM = std(PSTH_CS2All_CS_posi);
plotshaded([1:60],[mean(PSTH_CS2All_CS_posi)+SEM; mean(PSTH_CS2All_CS_posi); mean(PSTH_CS2All_CS_posi)-SEM]',[0.8 0 0],1)
SEM = std(PSTH_CS2All_CS_nega)/sqrt(length(PSTH_CS2All_CS_nega));
%SEM = std(PSTH_CS2All_CS_nega);
plotshaded([1:60],[mean(PSTH_CS2All_CS_nega)+SEM; mean(PSTH_CS2All_CS_nega); mean(PSTH_CS2All_CS_nega)-SEM]',[0 0 0.8],1)

%plot(mean(PSTH_CS2All_CS_posi),'r',LineWidth=3)
%plot(mean(PSTH_CS2All_CS_nega),'b',LineWidth=3)
plot([30 30],[-4 4],'--')
plot([0 60],[0 0],'--')
title('CS significant')
ylim([-4 4])
ylabel('z-scored(FR)')

subplot(2,2,2)
for ii=1:length(Opto_posi)
    PSTH_opto20Hz2All_Opto_posi = [PSTH_opto20Hz2All_Opto_posi; mean(PSTH_opto20Hz2Allz{Opto_posi(ii)})]; %z=scored
end

for ii=1:length(Opto_nega)
    PSTH_opto20Hz2All_Opto_nega = [PSTH_opto20Hz2All_Opto_nega; mean(PSTH_opto20Hz2Allz{Opto_nega(ii)})];
end


SEM = std(PSTH_opto20Hz2All_Opto_posi)/sqrt(length(PSTH_opto20Hz2All_Opto_posi));
%SEM = std(PSTH_opto20Hz2All_Opto_posi);
plotshaded([1:60],[mean(PSTH_opto20Hz2All_Opto_posi)+SEM; mean(PSTH_opto20Hz2All_Opto_posi); mean(PSTH_opto20Hz2All_Opto_posi)-SEM]',[0.8 0 0],1)
SEM = std(PSTH_opto20Hz2All_Opto_nega)/sqrt(length(PSTH_opto20Hz2All_Opto_nega));
%SEM = std(PSTH_opto20Hz2All_Opto_nega);
plotshaded([1:60],[mean(PSTH_opto20Hz2All_Opto_nega)+SEM; mean(PSTH_opto20Hz2All_Opto_nega); mean(PSTH_opto20Hz2All_Opto_nega)-SEM]',[0 0 0.8],1)

%plot(mean(PSTH_opto20Hz2All_Opto_posi),'r',LineWidth=3),hold on
%plot(mean(PSTH_opto20Hz2All_Opto_nega),'b',LineWidth=3)
plot([30 30],[-4 4],'--')
plot([0 60],[0 0],'--')
title('Opto significant')
ylim([-4 4])
ylabel('z-scored(FR)')

%
% {
PSTH_CS2All_tag=[];
PSTH_opto20Hz2All_tag=[];

subplot(2,2,3)
for ii=1:length(COI_tag)
    PSTH_CS2All_tag = [PSTH_CS2All_tag; mean(PSTH_CS2All{COI_tag(ii)})]; %firing rate
end

SEM = std(PSTH_CS2All_tag)/sqrt(length(PSTH_CS2All_tag));
plotshaded([1:60],[mean(PSTH_CS2All_tag)+SEM; mean(PSTH_CS2All_tag); mean(PSTH_CS2All_tag)-SEM]',[0.3 0.3 0.3],1)
plot([30 30],[0 6],'--')
title('Arc+ Tagged CS-Resp')
ylim([0 6])
ylabel('FR (Hz)')


subplot(2,2,4)
for ii=1:length(COI_tag)
    PSTH_opto20Hz2All_tag = [PSTH_opto20Hz2All_tag; mean(PSTH_opto20Hz2All{COI_tag(ii)})]; %firing rate
end

SEM = std(PSTH_opto20Hz2All_tag)/sqrt(length(PSTH_opto20Hz2All_tag));
plotshaded([1:60],[mean(PSTH_opto20Hz2All_tag)+SEM; mean(PSTH_opto20Hz2All_tag); mean(PSTH_opto20Hz2All_tag)-SEM]',[0.3 0.3 0.3],1)
plot([30 30],[0 6],'--')
title('Arc+ Tagged Opto-Resp')
ylim([0 6])
ylabel('FR (Hz)')
% }
%
PSTH_CS2All_All=[];
PSTH_opto20Hz2All_All=[];

subplot(2,2,1)
for ii=1:length(PSTH_CS2All)
    %plot(mean(PSTH_CS2All{ii}),Color=[0.5 0.5 0.5 0.2]),hold on
    %PSTH_CS2All_All = [PSTH_CS2All_All; mean(PSTH_CS2All{ii})];
    %plot(mean(PSTH_CS2Allz{ii}),Color=[0.5 0.5 0.5 0.2]),hold on
    PSTH_CS2All_All = [PSTH_CS2All_All; mean(PSTH_CS2Allz{ii})]; %z=scored
end

SEM = std(PSTH_CS2All_All)/sqrt(length(PSTH_CS2All_All));
%SEM = std(PSTH_CS2All_All);
plotshaded([1:60],[mean(PSTH_CS2All_All)+SEM; mean(PSTH_CS2All_All); mean(PSTH_CS2All_All)-SEM]',[0.5 0.5 0.5],1)
title('CS all')
ylim([-3 3])

subplot(2,2,2)

for ii=1:length(PSTH_opto20Hz2All)
    %plot(mean(PSTH_opto20Hz2All{ii}),Color=[0.5 0.5 0.5 0.2]),hold on
    %PSTH_opto20Hz2All_All = [PSTH_opto20Hz2All_All; mean(PSTH_opto20Hz2All{ii})];
    %plot(mean(PSTH_opto20Hz2Allz{ii}),Color=[0.5 0.5 0.5 0.2]),hold on
    PSTH_opto20Hz2All_All = [PSTH_opto20Hz2All_All; mean(PSTH_opto20Hz2Allz{ii})]; %z=scored
end

SEM = std(PSTH_opto20Hz2All_All)/sqrt(length(PSTH_opto20Hz2All_All));
%SEM = std(PSTH_opto20Hz2All_All);
plotshaded([1:60],[mean(PSTH_opto20Hz2All_All)+SEM; mean(PSTH_opto20Hz2All_All); mean(PSTH_opto20Hz2All_All)-SEM]',[0.5 0.5 0.5],1)
title('Opto all')
ylim([-3 3])

%%
baselineFR_tag=mean(PSTH_opto20Hz2All_tag(:,1:30),2);
OptoFR_tag=mean(PSTH_opto20Hz2All_tag(:,31:60),2);
%figure,
%errorbar(1,mean(baselineFR_tag),std(baselineFR_tag)/sqrt(length(baselineFR_tag))),hold on
%errorbar(2,mean(OptoFR_tag),std(OptoFR_tag)/sqrt(length(OptoFR_tag)))
%xlim([0.5 2.5])
%ylim([0 5])

%% Scatter regression?
%{
xx=mean(PSTH_CS2All_All(Opto_posi,31:40),2);
yy=mean(PSTH_opto20Hz2All_All(Opto_posi,31:40),2);

figure,subplot(1,3,1)
scatter(xx,yy,20,[0.9 0.5 0.5],"filled"),hold on
plot([0 0],[2 -2],'--k')
plot([-2 2],[0 0],'--k')

xx=mean(PSTH_CS2All_All(Opto_nega,31:40),2);
yy=mean(PSTH_opto20Hz2All_All(Opto_nega,31:40),2);

scatter(xx,yy,20,[0.5 0.5 0.9],"filled"),hold on
plot([0 0],[2 -2],'--k')
plot([-2 2],[0 0],'--k')

%
xx=mean(PSTH_CS2All_All(CS_posi,31:40),2);
yy=mean(PSTH_opto20Hz2All_All(CS_posi,31:40),2);

subplot(1,3,2)
scatter(xx,yy,20,[0.9 0.5 0.5],"filled"),hold on
plot([0 0],[2 -2],'--k')
plot([-2 2],[0 0],'--k')

xx=mean(PSTH_CS2All_All(CS_nega,31:40),2);
yy=mean(PSTH_opto20Hz2All_All(CS_nega,31:40),2);

scatter(xx,yy,20,[0.5 0.5 0.9],"filled"),hold on
plot([0 0],[2 -2],'--k')
plot([-2 2],[0 0],'--k')

%
xx=mean(PSTH_CS2All_All(:,31:40),2);
yy=mean(PSTH_opto20Hz2All_All(:,31:40),2);

subplot(1,3,3)
scatter(xx,yy,20,[0.5 0.5 0.5],"filled"),hold on
plot([0 0],[2 -2],'--k')
plot([-2 2],[0 0],'--k')

coefficients = polyfit(xx, yy, 1);

slope = coefficients(1);
intercept = coefficients(2);
plot(xx, slope * xx + intercept,'r')
%}
%% firing rate distribution
PSTH_CS2All_All=[];
PSTH_opto20Hz2All_All=[];

for ii=1:length(PSTH_CS2All)
    PSTH_CS2All_All = [PSTH_CS2All_All; mean(PSTH_CS2All{ii})];
    PSTH_opto20Hz2All_All = [PSTH_opto20Hz2All_All; mean(PSTH_opto20Hz2All{ii})];
end

CS_base=mean(PSTH_CS2All_All(:,1:30),2);
CS_sig=mean(PSTH_CS2All_All(:,31:60),2);
Opto_base=mean(PSTH_opto20Hz2All_All(:,1:30),2);
Opto_sig=mean(PSTH_opto20Hz2All_All(:,31:60),2);

figure,subplot(1,2,1)
edges = 10.^(-3:0.1:3);
scatter([1:60],histcounts(CS_base,edges),40,'k'),hold on
scatter([1:60],histcounts(CS_sig,edges),40,'r')
scatter([1:60],histcounts(Opto_sig,edges),40,'b')

%gaussian
f = @(x, params) params(1) * exp(-((x - params(2)).^2) / (2 * params(3)^2))
initialGuess = [100, 30, 30];

%fittedParams = lsqcurvefit(f, initialGuess, [1:60], histcounts(CS_base,edges));
fittedParams = fminsearch(@(params) norm(histcounts(CS_base,edges) - f([1:60], params))^2, initialGuess);
fittedCurve = f([1:60], fittedParams);
plot([1:60], fittedCurve, 'k');

fittedParams = fminsearch(@(params) norm(histcounts(CS_sig,edges) - f([1:60], params))^2, initialGuess);
fittedCurve = f([1:60], fittedParams);
plot([1:60], fittedCurve, 'r');

fittedParams = fminsearch(@(params) norm(histcounts(Opto_sig,edges) - f([1:60], params))^2, initialGuess);
fittedCurve = f([1:60], fittedParams);
plot([1:60], fittedCurve, 'b');

xticks([0 10 20 30 40 50 60])
xticklabels({'0.001','0.01','0.1','1','10','100','1000'})

subplot(1,2,2)
%Lorenz plot
[aa,bb]=sort(CS_base/max(CS_base));
plot(linspace(0,1,length(aa)),aa,'k'),hold on
[aa,bb]=sort(CS_sig/max(CS_sig));
plot(linspace(0,1,length(aa)),aa,'r'),hold on
[aa,bb]=sort(Opto_sig/max(Opto_sig));
plot(linspace(0,1,length(aa)),aa,'b'),hold on
grid on
set(gcf,'position',[500,500,1000,400])

%% spikewave t2peak
%{
for ii=1:1471
    wf_temp=DataUnitsAll.kilosort_wf{ii};
    [tr,tr_loc] = min(wf_temp);
    [peak,peak_loc] = max(wf_temp(tr_loc:end));
    t2p(ii)=peak_loc;
end

hist(t2p)

[hoge1,hoge2]=sort(t2p);

for ii=1:100
    figure,
    plot(DataUnitsAll.kilosort_wf{hoge2(ii)})
end
%}

%% animal analysis
for ii=1:10
    cellN(ii)=length(NselectAll{ii});
end
cellN=cumsum(cellN);

cellN_ID{1}=1:cellN(1);
for ii=2:10
    cellN_ID{ii}=cellN(ii-1)+1:cellN(ii);
end


% {
for ii=1:10
tmp=cellN_ID{ii};

tmp_tag=[];
tmp_exc=[];
tmp_inh=[];

for iii=1:length(tmp)
    if ismember(tmp(iii),COI_tag)
        tmp_tag=[tmp_tag tmp(iii)];
    elseif ismember(tmp(iii),COI_exc)
        tmp_exc=[tmp_exc tmp(iii)];
    elseif ismember(tmp(iii),COI_inh)
        tmp_inh=[tmp_inh tmp(iii)];
    end
end

cellN_tag{ii}=tmp_tag;
cellN_exc{ii}=tmp_exc;
cellN_inh{ii}=tmp_inh;
cellN_nc{ii}=setdiff(cellN_ID{ii},[cellN_tag{ii} cellN_exc{ii} cellN_inh{ii}]);

end
% }
PSTH_CS1All_All=[];
PSTH_opto20Hz1All_All=[];
for ii=1:length(PSTH_CS1All)
    PSTH_CS1All_All = [PSTH_CS1All_All; mean(PSTH_CS1All{ii})];
    PSTH_opto20Hz1All_All = [PSTH_opto20Hz1All_All; mean(PSTH_opto20Hz1All{ii})];
end


%% Mahalanobis Anal
CS=PSTH_CS1All_All(:,601:800);
Opto=PSTH_CS1All_All(:,401:600);
Base=PSTH_opto20Hz1All_All(:,401:600);

%partialdrop-tag
PT_Opto_tag=Opto;
PT_Base_tag=Base;
%PT_Opto_tag(COI_tag,:)=CS(COI_tag,:);
%PT_Base_tag(COI_tag,:)=CS(COI_tag,:);
%partialdrop-exc
PT_Opto_exc=Opto;
PT_Base_exc=Base;
%PT_Opto_exc(COI_exc,:)=CS(COI_exc,:);
%PT_Base_exc(COI_exc,:)=CS(COI_exc,:);
%partialdrop-inh
PT_Opto_inh=Opto;
PT_Base_inh=Base;
%PT_Opto_inh(COI_inh,:)=CS(COI_inh,:);
%PT_Base_inh(COI_inh,:)=CS(COI_inh,:);


for ii=1:10
d=mahalp(Opto(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); %mahalp(to-be-compared,mother)
d_CS_Opto(ii)=mean(d);
d=mahalp(Base(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_base(ii)=mean(d);
%{
d=mahalp(PT_Opto_tag(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_Opto_tag(ii)=mean(d);
%d=mahalp(PT_Base_tag(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d=mahalp(Base(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_base_tag(ii)=mean(d);

d=mahalp(PT_Opto_exc(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_Opto_exc(ii)=mean(d);
%d=mahalp(PT_Base_exc(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d=mahalp(Base(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_base_exc(ii)=mean(d);

d=mahalp(PT_Opto_inh(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_Opto_inh(ii)=mean(d);
%d=mahalp(PT_Base_inh(cellN_ID{ii},:)',CS(cellN_ID{ii},:)');
d=mahalp(Base(cellN_ID{ii},:)',CS(cellN_ID{ii},:)'); 
d_CS_base_inh(ii)=mean(d);
%}
% mean vector
mu_hat = mean(CS(cellN_ID{ii},:)');
sigma_hat = cov(CS(cellN_ID{ii},:)');

% Mahalanobis
mahalanobis_dist = mahalp(Opto(cellN_ID{ii},:)',CS(cellN_ID{ii},:)');

% contribution of each dim(cell)
dim_contrib = zeros(1, size(CS(cellN_ID{ii},:)', 2));
for i = 1:size(Base(cellN_ID{ii},:)', 2)
    data_temp = Base(cellN_ID{ii},:)';
    data_temp(:, i) = mu_hat(i); % replacing with mean
    mahalanobis_dist_temp = mahalp(data_temp, CS(cellN_ID{ii},:)');
    dim_contrib(i) = sum((mahalanobis_dist_temp - mahalanobis_dist).^2);
end

dim_contrib_normalized = dim_contrib / sum(dim_contrib);
contrbution_tag(ii)=sum(dim_contrib_normalized(cellN_tag{ii}-min(cellN_ID{ii})+1));
contrbution_exc(ii)=sum(dim_contrib_normalized(cellN_exc{ii}-min(cellN_ID{ii})+1));
contrbution_inh(ii)=sum(dim_contrib_normalized(cellN_inh{ii}-min(cellN_ID{ii})+1));

end

contrbution_nc = ones(1,10)-contrbution_tag-contrbution_exc-contrbution_inh;

%% contribution per cell
% { 
for ii=1:10
    if length(cellN_tag{ii})==0
        contrbution_tag_pcell(ii)=0;
    else
        contrbution_tag_pcell(ii)=contrbution_tag(ii)/length(cellN_tag{ii});
    end

    if length(cellN_exc{ii})==0
        contrbution_exc_pcell(ii)=0;
    else
        contrbution_exc_pcell(ii)=contrbution_exc(ii)/length(cellN_exc{ii});
    end

    if length(cellN_inh{ii})==0
        contrbution_inh_pcell(ii)=0;
    else
        contrbution_inh_pcell(ii)=contrbution_inh(ii)/length(cellN_inh{ii});
    end

    contrbution_nc_pcell(ii)=contrbution_nc(ii)/length(cellN_nc{ii});
end
% }
%% Drop Analysis Plot?
%{
figure,subplot(2,4,1)
plot([1 2],[d_CS_base; d_CS_Opto],'-o')
xlim([0.5 2.5])

subplot(2,4,5)
plot([1 2],[zeros(1,10); (d_CS_Opto-d_CS_base)./d_CS_base],'-',Color=[0.6 0.6 0.6]),hold on
plot([1 2],[0; mean((d_CS_Opto-d_CS_base)./d_CS_base)],'-',LineWidth=3,Color=[0.9 0.5 0.5])
plot([0 3], [0 0], '--')
xlim([0.5 2.5])

subplot(2,4,2)
plot([1 2],[d_CS_base_tag; d_CS_Opto_tag],'-o')
xlim([0.5 2.5])

subplot(2,4,6)
plot([1 2],[zeros(1,12); (d_CS_Opto_tag-d_CS_base_tag)./d_CS_base_tag],'-',Color=[0.6 0.6 0.6]),hold on
plot([1 2],[0; mean((d_CS_Opto_tag-d_CS_base_tag)./d_CS_base_tag)],'-',LineWidth=3,Color=[0.9 0.5 0.5])
plot([0 3], [0 0], '--')
xlim([0.5 2.5])

subplot(2,4,3)
plot([1 2],[d_CS_base_exc; d_CS_Opto_exc],'-o')
xlim([0.5 2.5])

subplot(2,4,7)
plot([1 2],[zeros(1,12); (d_CS_Opto_exc-d_CS_base_exc)./d_CS_base_exc],'-',Color=[0.6 0.6 0.6]),hold on
plot([1 2],[0; mean((d_CS_Opto_exc-d_CS_base_exc)./d_CS_base_exc)],'-',LineWidth=3,Color=[0.9 0.5 0.5])
plot([0 3], [0 0], '--')
xlim([0.5 2.5])

subplot(2,4,4)
plot([1 2],[d_CS_base_inh; d_CS_Opto_inh],'-o')
xlim([0.5 2.5])

subplot(2,4,8)
plot([1 2],[zeros(1,12); (d_CS_Opto_inh-d_CS_base_inh)./d_CS_base_inh],'-',Color=[0.6 0.6 0.6]),hold on
plot([1 2],[0; mean((d_CS_Opto_inh-d_CS_base_inh)./d_CS_base_inh)],'-',LineWidth=3,Color=[0.9 0.5 0.5])
plot([0 3], [0 0], '--')
xlim([0.5 2.5])
%}
%% Summary Plot
figure,
set(gcf,'position',[100,100,1200,400])
subplot(1,7,1)
%plot([1 2],[ones(1,12); (d_CS_Opto-d_CS_base)./d_CS_base],'-o')
plot([1 2],[zeros(1,10); (d_CS_Opto-d_CS_base)./d_CS_base],'-',Color=[0.6 0.6 0.6]),hold on
plot([1 2],[0; mean((d_CS_Opto-d_CS_base)./d_CS_base)],'-',LineWidth=3,Color=[0.9 0.5 0.5])
plot([0 3], [0 0], '--')
xlim([0.5 2.5])
% {
subplot(1,7,2)
bar([1,2,3,4],[mean(contrbution_tag) mean(contrbution_exc) mean(contrbution_inh) mean(contrbution_nc)]*100),hold on
errorbar([1,2,3,4],[mean(contrbution_tag) mean(contrbution_exc) mean(contrbution_inh) mean(contrbution_nc)]*100,[std(contrbution_tag) std(contrbution_exc) std(contrbution_inh) std(contrbution_nc)]*100/sqrt(12),'o')
scatter([1,2,3,4],[contrbution_tag; contrbution_exc; contrbution_inh; contrbution_nc]*100,[],[0.5 0.5 0.5],'filled')

subplot(1,7,3)
bar([1,2,3,4],[mean(contrbution_tag_pcell) mean(contrbution_exc_pcell) mean(contrbution_inh_pcell) mean(contrbution_nc_pcell)]*100),hold on
errorbar([1,2,3,4],[mean(contrbution_tag_pcell) mean(contrbution_exc_pcell) mean(contrbution_inh_pcell) mean(contrbution_nc_pcell)]*100,[std(contrbution_tag_pcell) std(contrbution_exc_pcell) std(contrbution_inh_pcell) std(contrbution_nc_pcell)]*100/sqrt(12),'o')
scatter([1,2,3,4],[contrbution_tag_pcell; contrbution_exc_pcell; contrbution_inh_pcell; contrbution_nc_pcell]*100,[],[0.5 0.5 0.5],'filled')
% }
% Freezing Anal

FreezingFCAll=[];
Ephys_FC{1}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc306_s01_fc\';
Ephys_FC{2}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc307_s01_fc\';
Ephys_FC{3}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc308_s01_fc\';
Ephys_FC{4}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc309_s01_fc\';
Ephys_FC{5}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc310_s01_fc\';
Ephys_FC{6}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc311_s01_fc\';
Ephys_FC{7}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc312_s01_fc\';
Ephys_FC{8}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc313_s01_fc\';
Ephys_FC{9}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc314_s01_fc\';
Ephys_FC{10}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc315_s01_fc\';

for ii=1:10
cd(Ephys_FC{ii})
load Freezing
FreezingFCAll=[FreezingFCAll; Freezing];
clear Freezing
end

FreezingFC=mean(FreezingFCAll(:,2:6),2);
dFreezingFC=mean(FreezingFCAll(:,2:6),2)-FreezingFCAll(:,1);
dFreezingFC_late=mean(FreezingFCAll(:,4:6),2)-FreezingFCAll(:,1);

FreezingAll=[];
Ephys_FC{1}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc306_s02_retrieval\';
Ephys_FC{2}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc307_s02_retrieval\';
Ephys_FC{3}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc308_s02_retrieval\';
Ephys_FC{4}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc309_s02_retrieval\';
Ephys_FC{5}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc310_s02_retrieval\';
Ephys_FC{6}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc311_s02_retrieval\';
Ephys_FC{7}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc312_s02_retrieval\';
Ephys_FC{8}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc313_s02_retrieval\';
Ephys_FC{9}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc314_s02_retrieval\';
Ephys_FC{10}='\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hagikent\IEGpaper\behavior\IEG_EphysBehavAnal\Arc315_s02_retrieval\';

for ii=1:10
cd(Ephys_FC{ii})
load Freezing
FreezingAll=[FreezingAll; Freezing];
clear Freezing
end

Freezing=mean(FreezingAll(:,2:6),2);
dFreezing=mean(FreezingAll(:,2:6),2)-FreezingAll(:,1);
dFreezing_early=mean(FreezingAll(:,2:3),2)-FreezingAll(:,1);

subplot(1,7,[4:5])
plot(1:6,FreezingFCAll','Color',[0.6 0.6 0.6]),hold on
errorbar(1:6,mean(FreezingFCAll),std(FreezingFCAll)/sqrt(12),LineWidth=3)

scatter(8,mean(FreezingAll(:,2:6),2)',50,[0.6 0.6 0.6],'filled'),hold on
errorbar(8,mean(mean(FreezingAll(:,2:6),2)),std(mean(FreezingAll(:,2:6),2))/sqrt(12),LineWidth=3)
xlim([0 9])
ylim([0 1])


%%
XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=FreezingFC;
YY=dFreezingFC_late;

subplot(1,7,[6:7])
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('FC freezing 1:5')


%% Correlation summary
figure,
set(gcf,'position',[100,100,2400,400])
%
subplot(1,6,1),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=FreezingFC;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('FC freezing 1:5')

%
subplot(1,6,2),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=dFreezingFC;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('FC freezing 1:5 - Base')

%
subplot(1,6,3),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=dFreezingFC_late;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('FC freezing 3:5 - Base')

%
subplot(1,6,4),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=Freezing;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('Retr. freezing 1:5')

%
subplot(1,6,5),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=dFreezing;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('Retr. freezing 1:5 - Base')

%
subplot(1,6,6),

XX=(d_CS_Opto-d_CS_base)./d_CS_base;
YY=dFreezing_early;
scatter(XX,YY,50,[0.5 0.5 0.5],'filled'),hold on

[p, R2, p_val,y_fit] = linear_fitKH(XX,YY');
R = corrcoef(XX,YY');
plot(XX, y_fit, 'b-',LineWidth=2)
ylim([-0.3 0.8])
title(['p = :' num2str(p_val) '; R = :' num2str(R(1,2))])
ylabel('Retr. freezing 1:2 - Base')


%% SVD

%for ii=1:10
ii=2;

CS_temp=CS(cellN_ID{ii},:);
Opto_temp=Opto(cellN_ID{ii},:);
Base_temp=Base(cellN_ID{ii},:);

[u,s,v_opto] = svd([Base_temp; Opto_temp]);
[u,s,rotation]=svd((Opto_temp-Base_temp)*v_opto(:,1:3));
v_opto(:,1:3) = v_opto(:,1:3)*rotation;

[u,s,v_CS] = svd([Base_temp; Opto_temp]);
[u,s,rotation]=svd((Opto_temp-Base_temp)*v_CS(:,1:3));
v_CS(:,1:3) = v_CS(:,1:3)*rotation;


figure(100),
plot(v_opto(:,2),v_opto(:,3)),hold on
plot(v_CS(:,2),v_CS(:,3))

%end

%%
%CS=PSTH_CS1All_All(:,601:800);
%Opto=PSTH_CS1All_All(:,401:600);
%Base=PSTH_opto20Hz1All_All(:,401:600);

CS=PSTH_CS1All_All(:,601:800);
Opto=PSTH_opto20Hz1All_All(:,601:800);
Base=PSTH_opto20Hz1All_All(:,401:600);

for ii=1:10

CS_temp=CS(cellN_ID{ii},:);
Opto_temp=Opto(cellN_ID{ii},:);
Base_temp=Base(cellN_ID{ii},:);

BaseAve_temp=mean(Base_temp,2);
CDCS = mean(CS_temp- BaseAve_temp,2);
Base_subtract = Base_temp - BaseAve_temp;

%CDCS length-square
CDCS_norm_squared = dot(CDCS, CDCS);

%Projecting Opto to CDCS
for jj=1:200
    projection_factor_opto(jj) = dot(Base_subtract(:,jj)', CDCS) / CDCS_norm_squared;
end

for jj=1:200
    projection_factor_opto(jj+200) = dot(Opto_temp(:,jj)', CDCS) / CDCS_norm_squared;
end

Opto_CDCS(:,ii)=projection_factor_opto;

end

%%
figure,
plot(Opto_CDCS,color = [0.7 0.7 0.7]),hold on
plot(mean(Opto_CDCS,2),color = [0.7 0.2 0.2],LineWidth=2)
ylim([-1,2])

