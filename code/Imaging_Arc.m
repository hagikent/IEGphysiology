%% Imaging_Arc
% code for total analysis

clear all
dir1='...\data\Imaging_Arc\1';
dir2='...\data\Imaging_Arc\2';
dir3='...\data\Imaging_Arc\3';
dir4='...\data\Imaging_Arc\4';
dir5='...\data\Imaging_Arc\5';
dir6='...\data\Imaging_Arc\6';
dir7='...\data\Imaging_Arc\7';

savedir=('...\data\Imaging_Arc');

cd(savedir)
load spontanstats.mat

Clims=[-0.0025 0.01];
LineColorGray=[0.9 0.9 0.9];
%%
PosiNn=[];
NegaNn=[];
cellN=[];

RespMat_Hab1_Posi=[];
RespMat_Hab2_Posi=[];
RespMat_FC_Posi=[];
RespMat_Ext1_Posi=[];
RespMat_Ext2_Posi=[];
RespMat_Test_Posi=[];

RespMat_Hab1_Nega=[];
RespMat_Hab2_Nega=[];
RespMat_FC_Nega=[];
RespMat_Ext1_Nega=[];
RespMat_Ext2_Nega=[];
RespMat_Test_Nega=[];

RespMat_Hab1_Un=[];
RespMat_Hab2_Un=[];
RespMat_FC_Un=[];
RespMat_Ext1_Un=[];
RespMat_Ext2_Un=[];
RespMat_Test_Un=[];

Hab1respAll_Posi=[];
Hab2respAll_Posi=[];
CSrespAll_Posi=[];
USrespAll_Posi=[];
Ext1respAll_Posi=[];
Ext2respAll_Posi=[];
TestrespAll_Posi=[];

Hab1respAll_Nega=[];
Hab2respAll_Nega=[];
CSrespAll_Nega=[];
USrespAll_Nega=[];
Ext1respAll_Nega=[];
Ext2respAll_Nega=[];
TestrespAll_Nega=[];

for ii=1:5
dir=eval(['dir',num2str(ii)]);
cd(dir)
cd('IJROItoTC_BLA')
load PosiN
load NegaN
cd(dir)
cd('TC_Anal')
load RespMat

RespMat_Hab1_Posi=cat(2,RespMat_Hab1_Posi,RespMat_Hab1(:,PosiN,:));
RespMat_Hab2_Posi=cat(2,RespMat_Hab2_Posi,RespMat_Hab2(:,PosiN,:));
RespMat_FC_Posi=cat(2,RespMat_FC_Posi,RespMat_FC(:,PosiN,:));
RespMat_Ext1_Posi=cat(2,RespMat_Ext1_Posi,RespMat_Ext1(:,PosiN,:));
RespMat_Ext2_Posi=cat(2,RespMat_Ext2_Posi,RespMat_Ext2(:,PosiN,:));
RespMat_Test_Posi=cat(2,RespMat_Test_Posi,RespMat_Test(:,PosiN,:));

RespMat_Hab1_Nega=cat(2,RespMat_Hab1_Nega,RespMat_Hab1(:,NegaN,:));
RespMat_Hab2_Nega=cat(2,RespMat_Hab2_Nega,RespMat_Hab2(:,NegaN,:));
RespMat_FC_Nega=cat(2,RespMat_FC_Nega,RespMat_FC(:,NegaN,:));
RespMat_Ext1_Nega=cat(2,RespMat_Ext1_Nega,RespMat_Ext1(:,NegaN,:));
RespMat_Ext2_Nega=cat(2,RespMat_Ext2_Nega,RespMat_Ext2(:,NegaN,:));
RespMat_Test_Nega=cat(2,RespMat_Test_Nega,RespMat_Test(:,NegaN,:));

RespMat_Hab1_Un=cat(2,RespMat_Hab1_Un,RespMat_Hab1(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Hab2_Un=cat(2,RespMat_Hab2_Un,RespMat_Hab2(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_FC_Un=cat(2,RespMat_FC_Un,RespMat_FC(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Ext1_Un=cat(2,RespMat_Ext1_Un,RespMat_Ext1(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Ext2_Un=cat(2,RespMat_Ext2_Un,RespMat_Ext2(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Test_Un=cat(2,RespMat_Test_Un,RespMat_Test(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));

PosiNn=[PosiNn size(PosiN,2)];
NegaNn=[NegaNn size(NegaN,2)];
PosiNall{ii}=PosiN;
NegaNall{ii}=NegaN;
cellN=[cellN size(RespMat_Hab1,2)];
load TCs

load resp
Hab1respAll_Posi=[Hab1respAll_Posi; Hab1resp(PosiN,:)];
Hab2respAll_Posi=[Hab2respAll_Posi; Hab2resp(PosiN,:)];
CSrespAll_Posi=[CSrespAll_Posi; CSresp(PosiN,:)];
USrespAll_Posi=[USrespAll_Posi; USresp(PosiN,:)];
Ext1respAll_Posi=[Ext1respAll_Posi; Ext1resp(PosiN,:)];
Ext2respAll_Posi=[Ext2respAll_Posi; Ext2resp(PosiN,:)];
TestrespAll_Posi=[TestrespAll_Posi; Testresp(PosiN,:)];

Hab1respAll_Nega=[Hab1respAll_Nega; Hab1resp(NegaN,:)];
Hab2respAll_Nega=[Hab2respAll_Nega; Hab2resp(NegaN,:)];
CSrespAll_Nega=[CSrespAll_Nega; CSresp(NegaN,:)];
USrespAll_Nega=[USrespAll_Nega; USresp(NegaN,:)];
Ext1respAll_Nega=[Ext1respAll_Nega; Ext1resp(NegaN,:)];
Ext2respAll_Nega=[Ext2respAll_Nega; Ext2resp(NegaN,:)];
TestrespAll_Nega=[TestrespAll_Nega; Testresp(NegaN,:)];

for jj=1:9
TCnormAll{ii,jj}=TCnorm{jj};
TCnorm_nonDTAll{ii,jj}=TCnorm_nonDT{jj};
TC_DTAll{ii,jj}=TC_DT{jj};
TC_nonDTAll{ii,jj}=TC_nonDT{jj};

temp=spontanstats{ii,jj};
for kk=1:length(PosiN)
spontanstats_Posi{ii,jj}{kk}=temp{PosiN(kk)};
end
for kk=1:length(NegaN)
spontanstats_Nega{ii,jj}{kk}=temp{NegaN(kk)};
end
clear temp
end

clear RespMat_Hab1 RespMat_Hab2 RespMat_FC RespMat_Ext1 RespMat_Ext2 RespMat_Test PosiN NegaN
clear TCnorm TCnorm_nonDT TC_DT TC_nonDT
clear Hab1resp Hab2resp CSresp USresp Ext1resp Ext2resp Testresp
end
%%
for ii=6:7
dir=eval(['dir',num2str(ii)]);
cd(dir)
load PosiN
load NegaN
cd('TC_Anal')
load RespMat

RespMat_Hab1_Posi=cat(2,RespMat_Hab1_Posi,RespMat_Hab1(:,PosiN,:));
RespMat_Hab2_Posi=cat(2,RespMat_Hab2_Posi,RespMat_Hab2(:,PosiN,:));
RespMat_FC_Posi=cat(2,RespMat_FC_Posi,RespMat_FC(:,PosiN,:));
RespMat_Ext1_Posi=cat(2,RespMat_Ext1_Posi,RespMat_Ext1(:,PosiN,:));
RespMat_Ext2_Posi=cat(2,RespMat_Ext2_Posi,RespMat_Ext2(:,PosiN,:));
RespMat_Test_Posi=cat(2,RespMat_Test_Posi,RespMat_Test(:,PosiN,:));

RespMat_Hab1_Nega=cat(2,RespMat_Hab1_Nega,RespMat_Hab1(:,NegaN,:));
RespMat_Hab2_Nega=cat(2,RespMat_Hab2_Nega,RespMat_Hab2(:,NegaN,:));
RespMat_FC_Nega=cat(2,RespMat_FC_Nega,RespMat_FC(:,NegaN,:));
RespMat_Ext1_Nega=cat(2,RespMat_Ext1_Nega,RespMat_Ext1(:,NegaN,:));
RespMat_Ext2_Nega=cat(2,RespMat_Ext2_Nega,RespMat_Ext2(:,NegaN,:));
RespMat_Test_Nega=cat(2,RespMat_Test_Nega,RespMat_Test(:,NegaN,:));

RespMat_Hab1_Un=cat(2,RespMat_Hab1_Un,RespMat_Hab1(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Hab2_Un=cat(2,RespMat_Hab2_Un,RespMat_Hab2(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_FC_Un=cat(2,RespMat_FC_Un,RespMat_FC(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Ext1_Un=cat(2,RespMat_Ext1_Un,RespMat_Ext1(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Ext2_Un=cat(2,RespMat_Ext2_Un,RespMat_Ext2(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));
RespMat_Test_Un=cat(2,RespMat_Test_Un,RespMat_Test(:,setdiff(1:size(RespMat_Hab1,2),union(PosiN,NegaN)),:));

PosiNn=[PosiNn size(PosiN,2)];
NegaNn=[NegaNn size(NegaN,2)];
PosiNall{ii}=PosiN;
NegaNall{ii}=NegaN;
cellN=[cellN size(RespMat_Hab1,2)];

load resp
Hab1respAll_Posi=[Hab1respAll_Posi; Hab1resp(PosiN,:)];
Hab2respAll_Posi=[Hab2respAll_Posi; Hab2resp(PosiN,:)];
CSrespAll_Posi=[CSrespAll_Posi; CSresp(PosiN,:)];
USrespAll_Posi=[USrespAll_Posi; USresp(PosiN,:)];
Ext1respAll_Posi=[Ext1respAll_Posi; Ext1resp(PosiN,:)];
Ext2respAll_Posi=[Ext2respAll_Posi; Ext2resp(PosiN,:)];
TestrespAll_Posi=[TestrespAll_Posi; Testresp(PosiN,:)];

Hab1respAll_Nega=[Hab1respAll_Nega; Hab1resp(NegaN,:)];
Hab2respAll_Nega=[Hab2respAll_Nega; Hab2resp(NegaN,:)];
CSrespAll_Nega=[CSrespAll_Nega; CSresp(NegaN,:)];
USrespAll_Nega=[USrespAll_Nega; USresp(NegaN,:)];
Ext1respAll_Nega=[Ext1respAll_Nega; Ext1resp(NegaN,:)];
Ext2respAll_Nega=[Ext2respAll_Nega; Ext2resp(NegaN,:)];
TestrespAll_Nega=[TestrespAll_Nega; Testresp(NegaN,:)];

load TCs
for jj=1:9
TCnormAll{ii,jj}=TCnorm{jj};
TCnorm_nonDTAll{ii,jj}=TCnorm_nonDT{jj};
TC_DTAll{ii,jj}=TC_DT{jj};
TC_nonDTAll{ii,jj}=TC_nonDT{jj};

temp=spontanstats{ii,jj};
for kk=1:length(PosiN)
spontanstats_Posi{ii,jj}{kk}=temp{PosiN(kk)};
end
for kk=1:length(NegaN)
spontanstats_Nega{ii,jj}{kk}=temp{NegaN(kk)};
end
clear temp
end

clear RespMat_Hab1 RespMat_Hab2 RespMat_FC RespMat_Ext1 RespMat_Ext2 RespMat_Test PosiN NegaN
clear TCnorm TCnorm_nonDT TC_DT TC_nonDT
clear Hab1resp Hab2resp CSresp USresp Ext1resp Ext2resp Testresp
end

%% baselining

RespMatBase_Hab1_Posi=repmat(mean(mean(RespMat_Hab1_Posi(1:50,:,1:5),3)),391,1);
RespMatBase_Hab2_Posi=repmat(mean(mean(RespMat_Hab2_Posi(1:50,:,1:5),3)),391,1);
RespMatBase_FC_Posi=repmat(mean(mean(RespMat_FC_Posi(1:50,:,1:5),3)),441,1);
RespMatBase_Ext1e_Posi=repmat(mean(mean(RespMat_Ext1_Posi(1:50,:,1:5),3)),391,1);
RespMatBase_Ext2e_Posi=repmat(mean(mean(RespMat_Ext2_Posi(1:50,:,1:5),3)),391,1);
RespMatBase_Ext1l_Posi=repmat(mean(mean(RespMat_Ext1_Posi(1:50,:,21:25),3)),391,1);
RespMatBase_Ext2l_Posi=repmat(mean(mean(RespMat_Ext2_Posi(1:50,:,21:25),3)),391,1);
RespMatBase_Test_Posi=repmat(mean(mean(RespMat_Test_Posi(1:50,:,1:5),3)),391,1);

RespMatBase_Hab1_Nega=repmat(mean(mean(RespMat_Hab1_Nega(1:50,:,1:5),3)),391,1);
RespMatBase_Hab2_Nega=repmat(mean(mean(RespMat_Hab2_Nega(1:50,:,1:5),3)),391,1);
RespMatBase_FC_Nega=repmat(mean(mean(RespMat_FC_Nega(1:50,:,1:5),3)),441,1);
RespMatBase_Ext1e_Nega=repmat(mean(mean(RespMat_Ext1_Nega(1:50,:,1:5),3)),391,1);
RespMatBase_Ext2e_Nega=repmat(mean(mean(RespMat_Ext2_Nega(1:50,:,1:5),3)),391,1);
RespMatBase_Ext1l_Nega=repmat(mean(mean(RespMat_Ext1_Nega(1:50,:,21:25),3)),391,1);
RespMatBase_Ext2l_Nega=repmat(mean(mean(RespMat_Ext2_Nega(1:50,:,21:25),3)),391,1);
RespMatBase_Test_Nega=repmat(mean(mean(RespMat_Test_Nega(1:50,:,1:5),3)),391,1);

RespMatBase_Hab1_Un=repmat(mean(mean(RespMat_Hab1_Un(1:50,:,1:5),3)),391,1);
RespMatBase_Hab2_Un=repmat(mean(mean(RespMat_Hab2_Un(1:50,:,1:5),3)),391,1);
RespMatBase_FC_Un=repmat(mean(mean(RespMat_FC_Un(1:50,:,1:5),3)),441,1);
RespMatBase_Ext1e_Un=repmat(mean(mean(RespMat_Ext1_Un(1:50,:,1:5),3)),391,1);
RespMatBase_Ext2e_Un=repmat(mean(mean(RespMat_Ext2_Un(1:50,:,1:5),3)),391,1);
RespMatBase_Ext1l_Un=repmat(mean(mean(RespMat_Ext1_Un(1:50,:,21:25),3)),391,1);
RespMatBase_Ext2l_Un=repmat(mean(mean(RespMat_Ext2_Un(1:50,:,21:25),3)),391,1);
RespMatBase_Test_Un=repmat(mean(mean(RespMat_Test_Un(1:50,:,1:5),3)),391,1);

%% self sorted
%Hab1

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Hab1_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Hab1_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Hab1_Un,2)
    temp=mean(RespMat_Hab1_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100),
f = figure(100);
f.Position = [100 100 2000 1500];
subplot(3,9,1),imagesc(mean(RespMat_Hab1_Posi(:,BBposi,:),3)' - RespMatBase_Hab1_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,10),imagesc(mean(RespMat_Hab1_Nega(:,BBnega,:),3)' - RespMatBase_Hab1_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,19),imagesc(mean(RespMat_Hab1_Un(:,BBun,:),3)' - RespMatBase_Hab1_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),
f = figure(101);
f.Position = [100 100 2000 200];
subplot(1,9,1)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%Hab2
for ii=1:size(RespMat_Hab2_Posi,2)
    temp=mean(RespMat_Hab2_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Hab2_Nega,2)
    temp=mean(RespMat_Hab2_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Hab2_Un,2)
    temp=mean(RespMat_Hab2_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100),
subplot(3,9,2),imagesc(mean(RespMat_Hab2_Posi(:,BBposi,:),3)'- RespMatBase_Hab2_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,11),imagesc(mean(RespMat_Hab2_Nega(:,BBnega,:),3)'- RespMatBase_Hab2_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,20),imagesc(mean(RespMat_Hab2_Un(:,BBun,:),3)'- RespMatBase_Hab2_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,2)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%FC-CS
for ii=1:size(RespMat_FC_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_FC_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_FC_Un,2)
    temp=mean(RespMat_FC_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100),
subplot(3,9,3),imagesc(mean(RespMat_FC_Posi(:,BBposi,:),3)'- RespMatBase_FC_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,12),imagesc(mean(RespMat_FC_Nega(:,BBnega,:),3)'- RespMatBase_FC_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,21),imagesc(mean(RespMat_FC_Un(:,BBun,:),3)'- RespMatBase_FC_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,3)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')


%FC-US
for ii=1:size(RespMat_FC_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_FC_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_FC_Un,2)
    temp=mean(RespMat_FC_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100),
subplot(3,9,4),imagesc(mean(RespMat_FC_Posi(:,BBposi,:),3)'- RespMatBase_FC_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,13),imagesc(mean(RespMat_FC_Nega(:,BBnega,:),3)'- RespMatBase_FC_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,22),imagesc(mean(RespMat_FC_Un(:,BBun,:),3)'- RespMatBase_FC_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,4)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%Ext1-5
for ii=1:size(RespMat_Ext1_Posi,2)
    temp=mean(RespMat_Ext1_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Ext1_Nega,2)
    temp=mean(RespMat_Ext1_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Ext1_Un,2)
    temp=mean(RespMat_Ext1_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp
figure(100),
subplot(3,9,5),imagesc(mean(RespMat_Ext1_Posi(:,BBposi,1:5),3)'- RespMatBase_Ext1e_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,14),imagesc(mean(RespMat_Ext1_Nega(:,BBnega,1:5),3)'- RespMatBase_Ext1e_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,23),imagesc(mean(RespMat_Ext1_Un(:,BBun,1:5),3)'- RespMatBase_Ext1e_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,5)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%Ext21-25
for ii=1:size(RespMat_Ext1_Posi,2)
    temp=mean(RespMat_Ext1_Posi(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Ext1_Nega,2)
    temp=mean(RespMat_Ext1_Nega(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Ext1_Un,2)
    temp=mean(RespMat_Ext1_Un(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100),
subplot(3,9,6),imagesc(mean(RespMat_Ext1_Posi(:,BBposi,21:25),3)'- RespMatBase_Ext1l_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,15),imagesc(mean(RespMat_Ext1_Nega(:,BBnega,21:25),3)'- RespMatBase_Ext1l_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,24),imagesc(mean(RespMat_Ext1_Un(:,BBun,21:25),3)'- RespMatBase_Ext1l_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,6)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%Ext2 1-5
for ii=1:size(RespMat_Ext2_Posi,2)
    temp=mean(RespMat_Ext2_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Ext2_Nega,2)
    temp=mean(RespMat_Ext2_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Ext2_Un,2)
    temp=mean(RespMat_Ext2_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100)
subplot(3,9,7),imagesc(mean(RespMat_Ext2_Posi(:,BBposi,1:5),3)'- RespMatBase_Ext2e_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,16),imagesc(mean(RespMat_Ext2_Nega(:,BBnega,1:5),3)'- RespMatBase_Ext2e_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,25),imagesc(mean(RespMat_Ext2_Un(:,BBun,1:5),3)'- RespMatBase_Ext2e_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,7)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')


%Ext2 21-25
for ii=1:size(RespMat_Ext2_Posi,2)
    temp=mean(RespMat_Ext2_Posi(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Ext2_Nega,2)
    temp=mean(RespMat_Ext2_Nega(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Ext2_Un,2)
    temp=mean(RespMat_Ext2_Un(:,ii,21:25),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100)
subplot(3,9,8),imagesc(mean(RespMat_Ext2_Posi(:,BBposi,21:25),3)'- RespMatBase_Ext2l_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,17),imagesc(mean(RespMat_Ext2_Nega(:,BBnega,21:25),3)'- RespMatBase_Ext2l_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,26),imagesc(mean(RespMat_Ext2_Un(:,BBun,21:25),3)'- RespMatBase_Ext2l_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,8)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

%ExtTest
for ii=1:size(RespMat_Test_Posi,2)
    temp=mean(RespMat_Test_Posi(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAposi,BBposi]=sort(resp,'descend');
clear resp

for ii=1:size(RespMat_Test_Nega,2)
    temp=mean(RespMat_Test_Nega(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAnega,BBnega]=sort(resp,'descend');
clear resp


for ii=1:size(RespMat_Test_Un,2)
    temp=mean(RespMat_Test_Un(:,ii,1:5),3);
    resp(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire US 4sec
end

[AAun,BBun]=sort(resp,'descend');
clear resp

figure(100)
subplot(3,9,9),imagesc(mean(RespMat_Test_Posi(:,BBposi,1:5),3)'- RespMatBase_Test_Posi(:,BBposi)',Clims),hold on, axis image
plot([50 50],[0 size(AAposi,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAposi,2)],'--','Color',[LineColorGray])
subplot(3,9,18),imagesc(mean(RespMat_Test_Nega(:,BBnega,1:5),3)'- RespMatBase_Test_Nega(:,BBnega)',Clims),hold on, axis image
plot([50 50],[0 size(AAnega,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAnega,2)],'--','Color',[LineColorGray])
subplot(3,9,27),imagesc(mean(RespMat_Test_Un(:,BBun,1:5),3)'- RespMatBase_Test_Un(:,BBun)',Clims),hold on, axis image
plot([50 50],[0 size(AAun,2)],'--','Color',LineColorGray)
plot([340 340],[0 size(AAun,2)],'--','Color',[LineColorGray])

figure(101),subplot(1,9,9)
plot(cumsum(histc(AAposi,[-0.02:0.0005:0.035]))/length(AAposi),'r'),hold on
plot(cumsum(histc(AAnega,[-0.02:0.0005:0.035]))/length(AAnega),'b')
plot(cumsum(histc(AAun,[-0.02:0.0005:0.035]))/length(AAun),'k')
plot([41 41],[0 1],'--','Color',[0.7 0.7 0.7])
xlim([0 100])
set(gca,'xtick',[0 41 80]); 
set(gca,'xticklabel',{'-2%','0','2%'});
xlabel('dF/F')
title('cumurative')
ylabel('fraction')

colormap(parula)


%% VS
for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    USrespPosi(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US + 3sec
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    USrespNega(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US + 3sec
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    FCCSrespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    FCCSrespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Ext1_Posi(:,ii,1:5),3);
    Ext1ArespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Ext1_Nega(:,ii,1:5),3);
    Ext1ArespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Test_Posi(:,ii,1:5),3);
    TestrespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Test_Nega(:,ii,1:5),3);
    TestrespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Hab2_Posi(:,ii,1:5),3);
    Hab2respPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Hab2_Nega(:,ii,1:5),3);
    Hab2respNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end



%%
figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,Ext1ArespNega,'b')
scatter(USrespPosi,Ext1ArespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('Ext1:1-5')

figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(Ext1ArespNega,TestrespNega,'b')
scatter(Ext1ArespPosi,TestrespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('Ext1:1-5'),ylabel('Test')

figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,FCCSrespNega,'b')
scatter(USrespPosi,FCCSrespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('CS(FC)')

figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,TestrespNega,'b')
scatter(USrespPosi,TestrespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('Test')

figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(Hab2respNega,Ext1ArespNega,'b')
scatter(Hab2respPosi,Ext1ArespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('Hab2'),ylabel('Ext1:1-5')

%% CS history
figure,
plot([1:5],[mean(Hab1respAll_Posi(:,1),1) mean(Hab1respAll_Posi(:,2),1) mean(Hab1respAll_Posi(:,3),1) mean(Hab1respAll_Posi(:,4),1) mean(Hab1respAll_Posi(:,5),1)]),hold on
plot([6:10],[mean(Hab2respAll_Posi(:,1),1) mean(Hab2respAll_Posi(:,2),1) mean(Hab2respAll_Posi(:,3),1) mean(Hab2respAll_Posi(:,4),1) mean(Hab2respAll_Posi(:,5),1)]),hold on
plot([11:15],[mean(CSrespAll_Posi(:,1),1) mean(CSrespAll_Posi(:,2),1) mean(CSrespAll_Posi(:,3),1) mean(CSrespAll_Posi(:,4),1) mean(CSrespAll_Posi(:,5),1)]),hold on
plot([16:20],[mean(mean(Ext1respAll_Posi(:,1:5),2)) mean(mean(Ext1respAll_Posi(:,6:10),2)) mean(mean(Ext1respAll_Posi(:,11:15),2))...
    mean(mean(Ext1respAll_Posi(:,16:20),2)) mean(mean(Ext1respAll_Posi(:,21:25),2))]),hold on
plot([21:25],[mean(mean(Ext2respAll_Posi(:,1:5),2)) mean(mean(Ext2respAll_Posi(:,6:10),2)) mean(mean(Ext2respAll_Posi(:,11:15),2))...
    mean(mean(Ext2respAll_Posi(:,16:20),2)) mean(mean(Ext2respAll_Posi(:,21:25),2))]),hold on

%% Spontan Quant
%{
figure(11),
%figure(12),
%figure(13),
%for bp=5
for bp=1:9
AUC1_Posi=[];
AUC2_Posi=[];
AUC3_Posi=[];
sda_Posi=[];
JVsp_Posi=[];

AUC1_Nega=[];
AUC2_Nega=[];
AUC3_Nega=[];
sda_Nega=[];
JVsp_Nega=[];

for ii=1:7
    for jj=1:size(spontanstats_Posi{ii, bp},2)
    temp(jj)=spontanstats_Posi{ii, bp}{1, jj}.AUC1;
    temp2(jj)=spontanstats_Posi{ii, bp}{1, jj}.AUC2;
    temp3(jj)=spontanstats_Posi{ii, bp}{1, jj}.AUC3;
    temp4(jj)=spontanstats_Posi{ii, bp}{1, jj}.sd_above;
    temp5(jj)=spontanstats_Posi{ii, bp}{1, jj}.JoVo_spikeN;    
    end
    
    AUC1_Posi=[AUC1_Posi temp]; 
    AUC2_Posi=[AUC2_Posi temp2];
    AUC3_Posi=[AUC3_Posi temp3];
    sda_Posi=[sda_Posi temp4];
    JVsp_Posi=[JVsp_Posi temp5];
    clear temp temp2 temp3 temp4 temp5
    
    for jj=1:size(spontanstats_Nega{ii, bp},2)
    temp(jj)=spontanstats_Nega{ii, bp}{1, jj}.AUC1;
    temp2(jj)=spontanstats_Nega{ii, bp}{1, jj}.AUC2;
    temp3(jj)=spontanstats_Nega{ii, bp}{1, jj}.AUC3;
    temp4(jj)=spontanstats_Nega{ii, bp}{1, jj}.sd_above;
    temp5(jj)=spontanstats_Nega{ii, bp}{1, jj}.JoVo_spikeN;    
    end
    AUC1_Nega=[AUC1_Nega temp]; 
    AUC2_Nega=[AUC2_Nega temp2];
    AUC3_Nega=[AUC3_Nega temp3];
    sda_Nega=[sda_Nega temp4];
    JVsp_Nega=[JVsp_Nega temp5];
    clear temp temp2 temp3 temp4 temp5
end
%
figure(11),
subplot(1,9,bp)
bar([1], mean(AUC1_Posi),'FaceColor',[0.8 0.4 0.4],'BarWidth',0.6),hold on
bar([2], mean(AUC1_Nega),'FaceColor',[0.4 0.4 0.8],'BarWidth',0.6)
e=errorbar([1 2], [mean(AUC1_Posi) mean(AUC1_Nega)], [std(AUC1_Posi)/sqrt(length(AUC1_Posi)) std(AUC1_Nega)/sqrt(length(AUC1_Nega))], '.');
e.Color = 'k';
plot(0.7,AUC1_Posi,'o','Color',[0.7 0.7 0.7 0.8])
plot(2.3,AUC1_Nega,'o','Color',[0.7 0.7 0.7 0.8])
xlim([0 3])
ylim([0 0.01])
pvalue_AUC1(bp)=ranksum(AUC1_Posi, AUC1_Nega);
    %{
    figure(12),
    subplot(1,9,bp)
    bar([1], mean(AUC2_Posi),'FaceColor',[0.8 0.6 0.6],'BarWidth',0.6),hold on
    bar([2], mean(AUC2_Nega),'FaceColor',[0.6 0.6 0.8],'BarWidth',0.6)
    e=errorbar([1 2], [mean(AUC2_Posi) mean(AUC2_Nega)], [std(AUC2_Posi)/sqrt(length(AUC2_Posi)) std(AUC2_Nega)/sqrt(length(AUC2_Nega))], '.');
    e.Color = 'k';
    plot(0.9,AUC2_Posi,'o','Color',[0.5 0.5 0.5])
    plot(2.1,AUC2_Nega,'o','Color',[0.5 0.5 0.5])
    xlim([0 3])
    ylim([0 0.01])
    pvalue_AUC2(bp)=ranksum(AUC2_Posi, AUC2_Nega);
    
    figure(13),
    subplot(1,9,bp)
    bar([1], mean(AUC3_Posi),'FaceColor',[0.8 0.6 0.6],'BarWidth',0.6),hold on
    bar([2], mean(AUC3_Nega),'FaceColor',[0.6 0.6 0.8],'BarWidth',0.6)
    e=errorbar([1 2], [mean(AUC3_Posi) mean(AUC3_Nega)], [std(AUC3_Posi)/sqrt(length(AUC3_Posi)) std(AUC3_Nega)/sqrt(length(AUC3_Nega))], '.');
    e.Color = 'k';
    plot(0.9,AUC3_Posi,'o','Color',[0.5 0.5 0.5])
    plot(2.1,AUC3_Nega,'o','Color',[0.5 0.5 0.5])
    xlim([0 3])
    ylim([0 0.01])
    pvalue_AUC3(bp)=ranksum(AUC3_Posi, AUC3_Nega);
    %}
end
%%
FCAUC_Posi=AUC1_Posi;
FCAUC_Nega=AUC1_Nega;
%HC1_Posi=AUC1_Posi;
%HC1_Nega=AUC1_Nega;

%HC2_Posi=AUC1_Posi;
%HC2_Nega=AUC1_Nega;

%HC3_Posi=AUC1_Posi;
%HC3_Nega=AUC1_Nega;


%%
for ii=1:length(AUC1_Nega)
    PN{1,ii}='Nega';
end
for ii=length(AUC1_Nega)+1:length(AUC1_Nega)+length(AUC1_Posi)
    PN{1,ii}='Posi';
end

figure,
scatterhist([AUC1_Nega AUC1_Posi],[USrespNega USrespPosi],'Group',PN','Color','br'),hold on
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6])
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6])
xlim([0 0.02]),ylim([-0.02 0.02])
xlabel('AUC(HC3)'),ylabel('US(FC)')
%
figure,
scatterhist([AUC1_Nega AUC1_Posi],[FCCSrespNega FCCSrespPosi],'Group',PN','Color','br'),hold on
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6])
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6])
xlim([0 0.02]),ylim([-0.02 0.02])
xlabel('AUC(HC3)'),ylabel('CS(FC)')

%
figure,
scatterhist([AUC1_Nega AUC1_Posi],[Ext1ArespNega Ext1ArespPosi],'Group',PN','Color','br'),hold on
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6])
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6])
xlim([0 0.02]),ylim([-0.02 0.02])
xlabel('AUC(HC3)'),ylabel('Ext1:1-5')
%%
figure,
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(AUC1_Nega,FCAUC_Nega,'b')
scatter(AUC1_Posi,FCAUC_Posi,'r')
xlim([0 0.02]),ylim([0 0.02])
xlabel('AUC(HC3)'),ylabel('FC_AUC')

%%
figure,
subplot(1,2,1)
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(HC3_Nega,HC1_Nega,'b')
scatter(HC3_Posi,HC1_Posi,'r')
xlabel('AUC(HC3)'),ylabel('AUC(HC1)')
xlim([0 0.015]),ylim([0 0.015]),
subplot(1,2,2)
plot([0 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(HC3_Nega,HC2_Nega,'b')
scatter(HC3_Posi,HC2_Posi,'r')
xlabel('AUC(HC3)'),ylabel('AUC(HC2)')
xlim([0 0.015]),ylim([0 0.015]),
%}


%% Spontan1 Corr
%{
for ll=1:9
for ii=1:7
temp_tc=TCnorm_nonDTAll{ii,ll}';
temp_PosiN=PosiNall{ii};
temp_NegaN=NegaNall{ii};

%order=[temp_PosiN temp_NegaN setdiff([1:size(temp_tc,1)],union(temp_PosiN,temp_NegaN))];
order=[temp_PosiN temp_NegaN];

%CorrMat=zeros(size(temp_tc,1),size(temp_tc,1));
CorrMat=zeros(length(order),length(order));
for jj=1:length(order)
    for kk=1:length(order)
        CF=corrcoef(temp_tc(order(jj),:),temp_tc(order(kk),:));
        CorrMat(jj,kk)=CF(1,2);
    end
end

CorrMatPosi=zeros(size(temp_PosiN,2),size(temp_PosiN,2));
for jj=1:size(temp_PosiN,2)
    for kk=1:size(temp_PosiN,2)
        CF=corrcoef(temp_tc(temp_PosiN(jj),:),temp_tc(temp_PosiN(kk),:));
        CorrMatPosi(jj,kk)=CF(1,2);
    end
end

CorrMatNega=zeros(size(temp_NegaN,2),size(temp_NegaN,2));
for jj=1:size(temp_NegaN,2)
    for kk=1:size(temp_NegaN,2)
        CF=corrcoef(temp_tc(temp_NegaN(jj),:),temp_tc(temp_NegaN(kk),:));
        CorrMatNega(jj,kk)=CF(1,2);
    end
end

CorrMatAll{ii}=CorrMat;
CorrMatPosiAll{ii}=CorrMatPosi;
CorrMatNegaAll{ii}=CorrMatNega;
%{
figure,
subplot(1,2,1),imagesc(CorrMatPosi,[-1 1])
subplot(1,2,2),imagesc(CorrMatNega,[-1 1])
%}
CorrDist(ii,:)=histcounts(CorrMat,[-1:0.1:1]);
CorrDistPosi(ii,:)=histcounts(CorrMatPosi,[-1:0.1:1]);
CorrDistNega(ii,:)=histcounts(CorrMatNega,[-1:0.1:1]);
%display(ii)

end

SumPosi=sum(CorrDistPosi,1);SumPosi(20)=0;
SumPosi=SumPosi/sum(SumPosi);
SumNega=sum(CorrDistNega,1);SumNega(20)=0;
SumNega=SumNega/sum(SumNega);
SumPosiAll{ll}=SumPosi;
SumNegaAll{ll}=SumNega;
clear CorrDist CorrDistPosi CorrDistNega
display(ll)
end
%%
%figure('position', [100 100 1200 600]),
figure,
for ll=1:9
subplot(2,9,ll),bar(SumPosiAll{ll},'FaceColor',[0.8 0.5 0.5]),xlim([0 20]),ylim([0 0.35]),
ax = gca; ax.XTick=(0:5:20); ax.XTickLabel = {'-1','-0.5','0','0.5', '1'};

subplot(2,9,ll+9),bar(SumNegaAll{ll},'FaceColor',[0.5 0.5 0.8]),xlim([0 20]),ylim([0 0.35]),
ax = gca; ax.XTick=(0:5:20); ax.XTickLabel = {'-1','-0.5','0','0.5', '1'};
end
%}

%% VS (201110transplanted)

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    USrespPosi(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US + 3sec
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    USrespNega(ii)=mean(temp(341:380))-mean(temp(1:49)); %entire US + 3sec
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:5),3);
    FCCSrespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:5),3);
    FCCSrespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Ext1_Posi(:,ii,1:5),3);
    Ext1ArespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Ext1_Nega(:,ii,1:5),3);
    Ext1ArespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Test_Posi(:,ii,1:5),3);
    TestrespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Test_Nega(:,ii,1:5),3);
    TestrespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Posi,2)
    temp=mean(RespMat_Hab2_Posi(:,ii,1:5),3);
    Hab2respPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_Hab2_Nega(:,ii,1:5),3);
    Hab2respNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end


for ii=1:size(RespMat_FC_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,1:2),3);
    FCCSErespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,1:2),3);
    FCCSErespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_FC_Posi,2)
    temp=mean(RespMat_FC_Posi(:,ii,3:5),3);
    FCCSLrespPosi(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

for ii=1:size(RespMat_Hab1_Nega,2)
    temp=mean(RespMat_FC_Nega(:,ii,3:5),3);
    FCCSLrespNega(ii)=mean(temp(50:340))-mean(temp(1:49)); %entire CS
end

%
f = figure;
f.Position = [1000 1000 1400 300];
subplot(1,4,1)
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,FCCSrespNega,'b')
scatter(USrespPosi,FCCSrespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('CS(FC)')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(USrespNega,FCCSrespNega);
plot(USrespNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(USrespPosi,FCCSrespPosi);
plot(USrespPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )

subplot(1,4,2)
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,FCCSLrespNega - FCCSErespNega,'b')
scatter(USrespPosi,FCCSLrespPosi - FCCSErespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('FC CS late - early')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(USrespNega,FCCSLrespNega - FCCSErespNega);
plot(USrespNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(USrespPosi,FCCSLrespPosi - FCCSErespPosi);
plot(USrespPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )

subplot(1,4,3)
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,Ext1ArespNega,'b')
scatter(USrespPosi,Ext1ArespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('Ext1:1-5')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(USrespNega,Ext1ArespNega);
plot(USrespNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(USrespPosi,Ext1ArespPosi);
plot(USrespPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )

subplot(1,4,4)
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(USrespNega,Ext1ArespNega - Hab2respNega,'b')
scatter(USrespPosi,Ext1ArespPosi - Hab2respPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('US'),ylabel('Ext1:1-5 - Hab2')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(USrespNega,Ext1ArespNega - Hab2respNega);
plot(USrespNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(USrespPosi,Ext1ArespPosi - Hab2respPosi);
plot(USrespPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )

%%
f = figure;
f.Position = [1000 1000 2500 500];
subplot(1,8,1:2),
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(Hab2respNega,Ext1ArespNega,'b')
scatter(Hab2respPosi,Ext1ArespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('Hab2'),ylabel('Ext1:1-5')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(Hab2respNega,Ext1ArespNega);
plot(Hab2respNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(Hab2respPosi,Ext1ArespPosi);
plot(Hab2respPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )


CSdiffPosi = Ext1ArespPosi-Hab2respPosi;
CSdiffNega = Ext1ArespNega-Hab2respNega;

subplot(1,8,3),
plot(0.9,CSdiffPosi,'o',color=[0.8 0.1 0.1]),hold on
plot(2.1,CSdiffNega,'o',color=[0.3 0.3 0.3]),hold on
plot([0 3],[0 0],'--')
xlim([0.5 2.5])
ylim([-0.01 0.02])

%length(find(CSdiffPosi>0))/length(CSdiffPosi)
%length(find(CSdiffNega>0))/length(CSdiffNega)

%errorbar(1, mean(CSdiffPosi), std(CSdiffPosi)/sqrt(length(CSdiffPosi)),'o',LineWidth=2)
%errorbar(2, mean(CSdiffNega), std(CSdiffNega)/sqrt(length(CSdiffNega)),'o',LineWidth=2)
boxchart(ones(length(CSdiffPosi),1)*1.2,CSdiffPosi, 'BoxWidth',0.3, MarkerStyle = 'none')
boxchart(ones(length(CSdiffNega),1)*1.8,CSdiffNega, 'BoxWidth',0.3, MarkerStyle = 'none')


subplot(1,8,5:6),
plot([0 0],[-0.02 0.02],'Color',[0.6 0.6 0.6]),hold on
plot([-0.02 0.02],[0 0],'Color',[0.6 0.6 0.6]),hold on
scatter(FCCSLrespNega - FCCSErespNega,Ext1ArespNega,'b')
scatter(FCCSLrespPosi - FCCSErespPosi,Ext1ArespPosi,'r')
xlim([-0.02 0.02]),ylim([-0.02 0.02])
xlabel('FC:CSlate- FC:CSearly'),ylabel('Ext1:1-5')

[pnega, R2nega, p_val,y_fit] = linear_fitKH(FCCSLrespNega - FCCSErespNega,Ext1ArespNega);
plot(FCCSLrespNega - FCCSErespNega, y_fit, 'b-',LineWidth=2)
[pposi, R2posi, p_val,y_fit] = linear_fitKH(FCCSLrespPosi - FCCSErespPosi,Ext1ArespPosi);
plot(FCCSLrespPosi - FCCSErespPosi, y_fit, 'r-',LineWidth=2)
legend(['R-squared:' num2str(R2nega) ' slope:' num2str(pnega(1))], ['Posi R-squared:' num2str(R2posi) 'Posi slope:' num2str(pposi(1))] )


%% 201202 added "significantly responsive"

base=1:50; %5sec base
%signal=51:240; %19sec tone sig
bin=10; %1sec bin

%US resp posi
TCall=mean(RespMat_FC_Posi(:,:,:),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_posiUS(ii) = TCranksum(TC, base, 341:380, bin);
    if mean(TC(base))<mean(TC(341:380))
        v_posiUS(ii)=1;
    else
        v_posiUS(ii)=-1;
    end
    vp_posiUS(ii)=p_posiUS(ii)*v_posiUS(ii);
end

%US resp nega
TCall=mean(RespMat_FC_Nega(:,:,:),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_negaUS(ii) = TCranksum(TC, base, 341:380, bin);
    if mean(TC(base))<mean(TC(341:380))
        v_negaUS(ii)=1;
    else
        v_negaUS(ii)=-1;
    end
    vp_negaUS(ii)=p_negaUS(ii)*v_negaUS(ii);
end

pUSexc=find(0<vp_posiUS & vp_posiUS<0.05);
nUSexc=find(0<vp_negaUS & vp_negaUS<0.05);

pUSrespSig = length(pUSexc)/size(RespMat_Hab1_Posi,2)
nUSrespSig = length(nUSexc)/size(RespMat_Hab1_Nega,2)

%% Kai2(231226)
total_A = 135;
positive_A = 29;
negative_A = total_A - positive_A;

total_B = 240;
positive_B = 37;
negative_B = total_B - positive_B;

observed_data = [positive_A, negative_A; positive_B, negative_B];

[p, Q] = chi2test(observed_data);
display(['shock responsive fraction during Retrieval; Fos+ vs Fos-: ' num2str(p)])
%% CS up/down (hab vs recall)

base=1:50; %5sec base
signal=51:240; %19sec tone sig
bin=10; %1sec bin

%Tone Resp Early (Hab2)
TCall=mean(RespMat_Hab2_Posi(:,:,:),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_posiE(ii) = TCranksum(TC, base, signal, bin);
    if mean(TC(base))<mean(TC(signal))
        v_posiE(ii)=1;
    else
        v_posiE(ii)=-1;
    end
    vp_posiE(ii)=p_posiE(ii)*v_posiE(ii);
end

TCall=mean(RespMat_Hab2_Nega(:,:,:),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_negaE(ii) = TCranksum(TC, base, signal, bin);
    if mean(TC(base))<mean(TC(signal))
        v_negaE(ii)=1;
    else
        v_negaE(ii)=-1;
    end
    vp_negaE(ii)=p_negaE(ii)*v_negaE(ii);
end

%Tone Resp Late (recall, Ext1,1-5)
TCall=mean(RespMat_Ext1_Posi(:,:,1:5),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_posiL(ii) = TCranksum(TC, base, signal, bin);
    if mean(TC(base))<mean(TC(signal))
        v_posiL(ii)=1;
    else
        v_posiL(ii)=-1;
    end
    vp_posiL(ii)=p_posiL(ii)*v_posiL(ii);
end

TCall=mean(RespMat_Ext1_Nega(:,:,1:5),3);
for ii=1:size(TCall,2)
    TC=TCall(:,ii);
    p_negaL(ii) = TCranksum(TC, base, signal, bin);
    if mean(TC(base))<mean(TC(signal))
        v_negaL(ii)=1;
    else
        v_negaL(ii)=-1;
    end
    vp_negaL(ii)=p_negaL(ii)*v_negaL(ii);
end

pEexc=find(0<vp_posiE & vp_posiE<0.05);
pLexc=find(0<vp_posiL & vp_posiL<0.05); %toneresponsive during recall
nEexc=find(0<vp_negaE & vp_negaE<0.05);
nLexc=find(0<vp_negaL & vp_negaL<0.05);  %toneresponsive during recall

pEinh=find(-0.05<vp_posiE & vp_posiE<0);
pLinh=find(-0.05<vp_posiL & vp_posiL<0);
nEinh=find(-0.05<vp_negaE & vp_negaE<0);
nLinh=find(-0.05<vp_negaL & vp_negaL<0);

pCSrecallSig = length(pLexc)/size(RespMat_Hab1_Posi,2)
nCSrecallSig = length(nLexc)/size(RespMat_Hab1_Nega,2)

%% Kai2(231226)
total_A = 135;
positive_A = 26;
negative_A = total_A - positive_A;

total_B = 240;
positive_B = 50;
negative_B = total_B - positive_B;

observed_data = [positive_A, negative_A; positive_B, negative_B];

[p, Q] = chi2test(observed_data);
display(['tone responsive fraction during Retrieval; Arc+ vs Arc-: ' num2str(p)])

%% CSup/CSdown
CSdiffPosi = Ext1ArespPosi-Hab2respPosi;
CSdiffNega = Ext1ArespNega-Hab2respNega;

figure,
plot(0.9,CSdiffPosi,'o',color=[0.8 0.1 0.1]),hold on
plot(2.1,CSdiffNega,'o',color=[0.3 0.3 0.3]),hold on
plot([0 3],[0 0],'--')
xlim([0.5 2.5])
ylim([-0.01 0.02])

%length(find(CSdiffPosi>0))/length(CSdiffPosi)
%length(find(CSdiffNega>0))/length(CSdiffNega)

%errorbar(1, mean(CSdiffPosi), std(CSdiffPosi)/sqrt(length(CSdiffPosi)),'o',LineWidth=2)
%errorbar(2, mean(CSdiffNega), std(CSdiffNega)/sqrt(length(CSdiffNega)),'o',LineWidth=2)
boxchart(ones(length(CSdiffPosi),1)*1.2,CSdiffPosi, MarkerStyle = 'none')
boxchart(ones(length(CSdiffNega),1)*1.8,CSdiffNega, MarkerStyle = 'none')

ranksum(CSdiffPosi,CSdiffNega)

%% 230610 ave traces

Hab1_Posi=mean(RespMat_Hab1_Posi,3)-RespMatBase_Hab1_Posi;
Hab1_Posi_ave=mean(Hab1_Posi,2);
Hab1_Posi_sem=std(Hab1_Posi,[],2)/sqrt(size(Hab1_Posi,2));
Hab1_Nega=mean(RespMat_Hab1_Nega,3)-RespMatBase_Hab1_Nega;
Hab1_Nega_ave=mean(Hab1_Nega,2);
Hab1_Nega_sem=std(Hab1_Nega,[],2)/sqrt(size(Hab1_Nega,2));

Hab2_Posi=mean(RespMat_Hab2_Posi,3)-RespMatBase_Hab2_Posi;
Hab2_Posi_ave=mean(Hab2_Posi,2);
Hab2_Posi_sem=std(Hab2_Posi,[],2)/sqrt(size(Hab2_Posi,2));
Hab2_Nega=mean(RespMat_Hab2_Nega,3)-RespMatBase_Hab2_Nega;
Hab2_Nega_ave=mean(Hab2_Nega,2);
Hab2_Nega_sem=std(Hab2_Nega,[],2)/sqrt(size(Hab2_Nega,2));

FC_Posi=mean(RespMat_FC_Posi,3)-RespMatBase_FC_Posi;
FC_Posi_ave=mean(FC_Posi,2);
FC_Posi_sem=std(FC_Posi,[],2)/sqrt(size(FC_Posi,2));
FC_Nega=mean(RespMat_FC_Nega,3)-RespMatBase_FC_Nega;
FC_Nega_ave=mean(FC_Nega,2);
FC_Nega_sem=std(FC_Nega,[],2)/sqrt(size(FC_Nega,2));

Ext1e_Posi=mean(RespMat_Ext1_Posi(:,:,1:5),3)-RespMatBase_Ext1e_Posi;
Ext1e_Posi_ave=mean(Ext1e_Posi,2);
Ext1e_Posi_sem=std(Ext1e_Posi,[],2)/sqrt(size(Ext1e_Posi,2));
Ext1e_Nega=mean(RespMat_Ext1_Nega(:,:,1:5),3)-RespMatBase_Ext1e_Nega;
Ext1e_Nega_ave=mean(Ext1e_Nega,2);
Ext1e_Nega_sem=std(Ext1e_Nega,[],2)/sqrt(size(Ext1e_Nega,2));

Ext1l_Posi=mean(RespMat_Ext1_Posi(:,:,21:25),3)-RespMatBase_Ext1l_Posi;
Ext1l_Posi_ave=mean(Ext1l_Posi,2);
Ext1l_Posi_sem=std(Ext1l_Posi,[],2)/sqrt(size(Ext1l_Posi,2));
Ext1l_Nega=mean(RespMat_Ext1_Nega(:,:,21:25),3)-RespMatBase_Ext1l_Nega;
Ext1l_Nega_ave=mean(Ext1l_Nega,2);
Ext1l_Nega_sem=std(Ext1l_Nega,[],2)/sqrt(size(Ext1l_Nega,2));

Ext2e_Posi=mean(RespMat_Ext2_Posi(:,:,1:5),3)-RespMatBase_Ext2e_Posi;
Ext2e_Posi_ave=mean(Ext2e_Posi,2);
Ext2e_Posi_sem=std(Ext2e_Posi,[],2)/sqrt(size(Ext2e_Posi,2));
Ext2e_Nega=mean(RespMat_Ext2_Nega(:,:,1:5),3)-RespMatBase_Ext2e_Nega;
Ext2e_Nega_ave=mean(Ext2e_Nega,2);
Ext2e_Nega_sem=std(Ext2e_Nega,[],2)/sqrt(size(Ext2e_Nega,2));

Ext2l_Posi=mean(RespMat_Ext2_Posi(:,:,21:25),3)-RespMatBase_Ext2l_Posi;
Ext2l_Posi_ave=mean(Ext2l_Posi,2);
Ext2l_Posi_sem=std(Ext2l_Posi,[],2)/sqrt(size(Ext2l_Posi,2));
Ext2l_Nega=mean(RespMat_Ext2_Nega(:,:,21:25),3)-RespMatBase_Ext2l_Nega;
Ext2l_Nega_ave=mean(Ext2l_Nega,2);
Ext2l_Nega_sem=std(Ext2l_Nega,[],2)/sqrt(size(Ext2l_Nega,2));

Test_Posi=mean(RespMat_Test_Posi,3)-RespMatBase_Test_Posi;
Test_Posi_ave=mean(Test_Posi,2);
Test_Posi_sem=std(Test_Posi,[],2)/sqrt(size(Test_Posi,2));
Test_Nega=mean(RespMat_Test_Nega,3)-RespMatBase_Test_Nega;
Test_Nega_ave=mean(Test_Nega,2);
Test_Nega_sem=std(Test_Nega,[],2)/sqrt(size(Test_Nega,2));

%AUC (230624)
Hab1_Posi_AUC=mean(Hab1_Posi(51:340,:));
Hab2_Posi_AUC=mean(Hab2_Posi(51:340,:));
FCCS_Posi_AUC=mean(FC_Posi(51:340,:));
FCUS_Posi_AUC=mean(FC_Posi(340:380,:));
Ext1e_Posi_AUC=mean(Ext1e_Posi(51:340,:));
Ext1l_Posi_AUC=mean(Ext1l_Posi(51:340,:));
Ext2e_Posi_AUC=mean(Ext2e_Posi(51:340,:));
Ext2l_Posi_AUC=mean(Ext2l_Posi(51:340,:));
Test_Posi_AUC=mean(Test_Posi(51:340,:));

Hab1_Nega_AUC=mean(Hab1_Nega(51:340,:));
Hab2_Nega_AUC=mean(Hab2_Nega(51:340,:));
FCCS_Nega_AUC=mean(FC_Nega(51:340,:));
FCUS_Nega_AUC=mean(FC_Nega(340:380,:));
Ext1e_Nega_AUC=mean(Ext1e_Nega(51:340,:));
Ext1l_Nega_AUC=mean(Ext1l_Nega(51:340,:));
Ext2e_Nega_AUC=mean(Ext2e_Nega(51:340,:));
Ext2l_Nega_AUC=mean(Ext2l_Nega(51:340,:));
Test_Nega_AUC=mean(Test_Nega(51:340,:));

%%
YYlim=[-0.001 0.015];

f = figure(123);
f.Position = [100 100 2400 1500];

subplot(3,8,1)
plotshaded([1:size(Hab1_Posi,1)],[Hab1_Posi_ave+Hab1_Posi_sem Hab1_Posi_ave Hab1_Posi_ave-Hab1_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab1_Nega,1)],[Hab1_Nega_ave+Hab1_Nega_sem Hab1_Nega_ave Hab1_Nega_ave-Hab1_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab1')

subplot(3,8,2)
plotshaded([1:size(Hab2_Posi,1)],[Hab2_Posi_ave+Hab2_Posi_sem Hab2_Posi_ave Hab2_Posi_ave-Hab2_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab2_Nega,1)],[Hab2_Nega_ave+Hab2_Nega_sem Hab2_Nega_ave Hab2_Nega_ave-Hab2_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab2')

subplot(3,8,3)
plotshaded([1:size(FC_Posi,1)],[FC_Posi_ave+FC_Posi_sem FC_Posi_ave FC_Posi_ave-FC_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(FC_Nega,1)],[FC_Nega_ave+FC_Nega_sem FC_Nega_ave FC_Nega_ave-FC_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('FC')

subplot(3,8,4)
plotshaded([1:size(Ext1e_Posi,1)],[Ext1e_Posi_ave+Ext1e_Posi_sem Ext1e_Posi_ave Ext1e_Posi_ave-Ext1e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1e_Nega,1)],[Ext1e_Nega_ave+Ext1e_Nega_sem Ext1e_Nega_ave Ext1e_Nega_ave-Ext1e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial1-5')

subplot(3,8,5)
plotshaded([1:size(Ext1l_Posi,1)],[Ext1l_Posi_ave+Ext1l_Posi_sem Ext1l_Posi_ave Ext1l_Posi_ave-Ext1l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1l_Nega,1)],[Ext1l_Nega_ave+Ext1l_Nega_sem Ext1l_Nega_ave Ext1l_Nega_ave-Ext1l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial21-25')

subplot(3,8,6)
plotshaded([1:size(Ext2e_Posi,1)],[Ext2e_Posi_ave+Ext2e_Posi_sem Ext2e_Posi_ave Ext2e_Posi_ave-Ext2e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2e_Nega,1)],[Ext2e_Nega_ave+Ext2e_Nega_sem Ext2e_Nega_ave Ext2e_Nega_ave-Ext2e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial1-5')

subplot(3,8,7)
plotshaded([1:size(Ext2l_Posi,1)],[Ext2l_Posi_ave+Ext2l_Posi_sem Ext2l_Posi_ave Ext2l_Posi_ave-Ext2l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2l_Nega,1)],[Ext2l_Nega_ave+Ext2l_Nega_sem Ext2l_Nega_ave Ext2l_Nega_ave-Ext2l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial21-25')

subplot(3,8,8)
plotshaded([1:size(Test_Posi,1)],[Test_Posi_ave+Ext2l_Posi_sem Test_Posi_ave Test_Posi_ave-Test_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Test_Nega,1)],[Test_Nega_ave+Test_Nega_sem Test_Nega_ave Test_Nega_ave-Test_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('ExtTest')

%AUC (230624)
f = figure(124);
f.Position = [100 100 2400 1500];
subplot(3,9,1)
scatter(0.9,Hab1_Posi_AUC,'k'),hold on,scatter(2.1,Hab1_Nega_AUC,'k')
errorbar(1, mean(Hab1_Posi_AUC), std(Hab1_Posi_AUC)/sqrt(length(Hab1_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab1_Nega_AUC), std(Hab1_Nega_AUC)/sqrt(length(Hab1_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab1  p=' num2str(ranksum(Hab1_Posi_AUC, Hab1_Nega_AUC))])

subplot(3,9,2)
scatter(0.9,Hab2_Posi_AUC,'k'),hold on,scatter(2.1,Hab2_Nega_AUC,'k')
errorbar(1, mean(Hab2_Posi_AUC), std(Hab2_Posi_AUC)/sqrt(length(Hab2_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab2_Nega_AUC), std(Hab2_Nega_AUC)/sqrt(length(Hab2_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab2  p=' num2str(ranksum(Hab2_Posi_AUC, Hab2_Nega_AUC))])

subplot(3,9,3)
scatter(0.9,FCCS_Posi_AUC,'k'),hold on,scatter(2.1,FCCS_Nega_AUC,'k')
errorbar(1, mean(FCCS_Posi_AUC), std(FCCS_Posi_AUC)/sqrt(length(FCCS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCCS_Nega_AUC), std(FCCS_Nega_AUC)/sqrt(length(FCCS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC-CS  p=' num2str(ranksum(FCCS_Posi_AUC, FCCS_Nega_AUC))])

subplot(3,9,4)
scatter(0.9,FCUS_Posi_AUC,'k'),hold on,scatter(2.1,FCUS_Nega_AUC,'k')
errorbar(1, mean(FCUS_Posi_AUC), std(FCUS_Posi_AUC)/sqrt(length(FCUS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCUS_Nega_AUC), std(FCUS_Nega_AUC)/sqrt(length(FCUS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC-US  p=' num2str(ranksum(FCUS_Posi_AUC, FCUS_Nega_AUC))])

subplot(3,9,5)
scatter(0.9,Ext1e_Posi_AUC,'k'),hold on,scatter(2.1,Ext1e_Nega_AUC,'k')
errorbar(1, mean(Ext1e_Posi_AUC), std(Ext1e_Posi_AUC)/sqrt(length(Ext1e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1e_Nega_AUC), std(Ext1e_Nega_AUC)/sqrt(length(Ext1e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1 1-5  p=' num2str(ranksum(Ext1e_Posi_AUC, Ext1e_Nega_AUC))])

subplot(3,9,6)
scatter(0.9,Ext1l_Posi_AUC,'k'),hold on,scatter(2.1,Ext1l_Nega_AUC,'k')
errorbar(1, mean(Ext1l_Posi_AUC), std(Ext1l_Posi_AUC)/sqrt(length(Ext1l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1l_Nega_AUC), std(Ext1l_Nega_AUC)/sqrt(length(Ext1l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1_1-5  p=' num2str(ranksum(Ext1l_Posi_AUC, Ext1l_Nega_AUC))])

subplot(3,9,7)
scatter(0.9,Ext2e_Posi_AUC,'k'),hold on,scatter(2.1,Ext2e_Nega_AUC,'k')
errorbar(1, mean(Ext2e_Posi_AUC), std(Ext2e_Posi_AUC)/sqrt(length(Ext2e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2e_Nega_AUC), std(Ext2e_Nega_AUC)/sqrt(length(Ext2e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2_1-5  p=' num2str(ranksum(Ext2e_Posi_AUC, Ext2e_Nega_AUC))])

subplot(3,9,8)
scatter(0.9,Ext2l_Posi_AUC,'k'),hold on,scatter(2.1,Ext2l_Nega_AUC,'k')
errorbar(1, mean(Ext2l_Posi_AUC), std(Ext2l_Posi_AUC)/sqrt(length(Ext2l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2l_Nega_AUC), std(Ext2l_Nega_AUC)/sqrt(length(Ext2l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2 21-25  p=' num2str(ranksum(Ext2l_Posi_AUC, Ext2l_Nega_AUC))])

subplot(3,9,9)
scatter(0.9,Test_Posi_AUC,'k'),hold on,scatter(2.1,Test_Nega_AUC,'k')
errorbar(1, mean(Test_Posi_AUC), std(Test_Posi_AUC)/sqrt(length(Test_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Test_Nega_AUC), std(Test_Nega_AUC)/sqrt(length(Test_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Test  p=' num2str(ranksum(Test_Posi_AUC, Test_Nega_AUC))])


%% 230610 ave traces (CS responsive during recall)

Hab1_Posi=mean(RespMat_Hab1_Posi,3)-RespMatBase_Hab1_Posi;
Hab1_Posi=Hab1_Posi(:,pLexc);
Hab1_Posi_ave=mean(Hab1_Posi,2);
Hab1_Posi_sem=std(Hab1_Posi,[],2)/sqrt(size(Hab1_Posi,2));
Hab1_Nega=mean(RespMat_Hab1_Nega,3)-RespMatBase_Hab1_Nega;
Hab1_Nega=Hab1_Nega(:,nLexc);
Hab1_Nega_ave=mean(Hab1_Nega,2);
Hab1_Nega_sem=std(Hab1_Nega,[],2)/sqrt(size(Hab1_Nega,2));

Hab2_Posi=mean(RespMat_Hab2_Posi,3)-RespMatBase_Hab2_Posi;
Hab2_Posi=Hab2_Posi(:,pLexc);
Hab2_Posi_ave=mean(Hab2_Posi,2);
Hab2_Posi_sem=std(Hab2_Posi,[],2)/sqrt(size(Hab2_Posi,2));
Hab2_Nega=mean(RespMat_Hab2_Nega,3)-RespMatBase_Hab2_Nega;
Hab2_Nega=Hab2_Nega(:,nLexc);
Hab2_Nega_ave=mean(Hab2_Nega,2);
Hab2_Nega_sem=std(Hab2_Nega,[],2)/sqrt(size(Hab2_Nega,2));

FC_Posi=mean(RespMat_FC_Posi,3)-RespMatBase_FC_Posi;
FC_Posi=FC_Posi(:,pLexc);
FC_Posi_ave=mean(FC_Posi,2);
FC_Posi_sem=std(FC_Posi,[],2)/sqrt(size(FC_Posi,2));
FC_Nega=mean(RespMat_FC_Nega,3)-RespMatBase_FC_Nega;
FC_Nega=FC_Nega(:,nLexc);
FC_Nega_ave=mean(FC_Nega,2);
FC_Nega_sem=std(FC_Nega,[],2)/sqrt(size(FC_Nega,2));

Ext1e_Posi=mean(RespMat_Ext1_Posi(:,:,1:5),3)-RespMatBase_Ext1e_Posi;
Ext1e_Posi=Ext1e_Posi(:,pLexc);
Ext1e_Posi_ave=mean(Ext1e_Posi,2);
Ext1e_Posi_sem=std(Ext1e_Posi,[],2)/sqrt(size(Ext1e_Posi,2));
Ext1e_Nega=mean(RespMat_Ext1_Nega(:,:,1:5),3)-RespMatBase_Ext1e_Nega;
Ext1e_Nega=Ext1e_Nega(:,nLexc);
Ext1e_Nega_ave=mean(Ext1e_Nega,2);
Ext1e_Nega_sem=std(Ext1e_Nega,[],2)/sqrt(size(Ext1e_Nega,2));

Ext1l_Posi=mean(RespMat_Ext1_Posi(:,:,21:25),3)-RespMatBase_Ext1l_Posi;
Ext1l_Posi=Ext1l_Posi(:,pLexc);
Ext1l_Posi_ave=mean(Ext1l_Posi,2);
Ext1l_Posi_sem=std(Ext1l_Posi,[],2)/sqrt(size(Ext1l_Posi,2));
Ext1l_Nega=mean(RespMat_Ext1_Nega(:,:,21:25),3)-RespMatBase_Ext1l_Nega;
Ext1l_Nega=Ext1l_Nega(:,nLexc);
Ext1l_Nega_ave=mean(Ext1l_Nega,2);
Ext1l_Nega_sem=std(Ext1l_Nega,[],2)/sqrt(size(Ext1l_Nega,2));

Ext2e_Posi=mean(RespMat_Ext2_Posi(:,:,1:5),3)-RespMatBase_Ext2e_Posi;
Ext2e_Posi=Ext2e_Posi(:,pLexc);
Ext2e_Posi_ave=mean(Ext2e_Posi,2);
Ext2e_Posi_sem=std(Ext2e_Posi,[],2)/sqrt(size(Ext2e_Posi,2));
Ext2e_Nega=mean(RespMat_Ext2_Nega(:,:,1:5),3)-RespMatBase_Ext2e_Nega;
Ext2e_Nega=Ext2e_Nega(:,nLexc);
Ext2e_Nega_ave=mean(Ext2e_Nega,2);
Ext2e_Nega_sem=std(Ext2e_Nega,[],2)/sqrt(size(Ext2e_Nega,2));

Ext2l_Posi=mean(RespMat_Ext2_Posi(:,:,21:25),3)-RespMatBase_Ext2l_Posi;
Ext2l_Posi=Ext2l_Posi(:,pLexc);
Ext2l_Posi_ave=mean(Ext2l_Posi,2);
Ext2l_Posi_sem=std(Ext2l_Posi,[],2)/sqrt(size(Ext2l_Posi,2));
Ext2l_Nega=mean(RespMat_Ext2_Nega(:,:,21:25),3)-RespMatBase_Ext2l_Nega;
Ext2l_Nega=Ext2l_Nega(:,nLexc);
Ext2l_Nega_ave=mean(Ext2l_Nega,2);
Ext2l_Nega_sem=std(Ext2l_Nega,[],2)/sqrt(size(Ext2l_Nega,2));

Test_Posi=mean(RespMat_Test_Posi,3)-RespMatBase_Test_Posi;
Test_Posi=Test_Posi(:,pLexc);
Test_Posi_ave=mean(Test_Posi,2);
Test_Posi_sem=std(Test_Posi,[],2)/sqrt(size(Test_Posi,2));
Test_Nega=mean(RespMat_Test_Nega,3)-RespMatBase_Test_Nega;
Test_Nega=Test_Nega(:,nLexc);
Test_Nega_ave=mean(Test_Nega,2);
Test_Nega_sem=std(Test_Nega,[],2)/sqrt(size(Test_Nega,2));

%AUC (230624)
Hab1_Posi_AUC=mean(Hab1_Posi(51:340,:));
Hab2_Posi_AUC=mean(Hab2_Posi(51:340,:));
FCCS_Posi_AUC=mean(FC_Posi(51:340,:));
FCUS_Posi_AUC=mean(FC_Posi(340:380,:));
Ext1e_Posi_AUC=mean(Ext1e_Posi(51:340,:));
Ext1l_Posi_AUC=mean(Ext1l_Posi(51:340,:));
Ext2e_Posi_AUC=mean(Ext2e_Posi(51:340,:));
Ext2l_Posi_AUC=mean(Ext2l_Posi(51:340,:));
Test_Posi_AUC=mean(Test_Posi(51:340,:));

Hab1_Nega_AUC=mean(Hab1_Nega(51:340,:));
Hab2_Nega_AUC=mean(Hab2_Nega(51:340,:));
FCCS_Nega_AUC=mean(FC_Nega(51:340,:));
FCUS_Nega_AUC=mean(FC_Nega(340:380,:));
Ext1e_Nega_AUC=mean(Ext1e_Nega(51:340,:));
Ext1l_Nega_AUC=mean(Ext1l_Nega(51:340,:));
Ext2e_Nega_AUC=mean(Ext2e_Nega(51:340,:));
Ext2l_Nega_AUC=mean(Ext2l_Nega(51:340,:));
Test_Nega_AUC=mean(Test_Nega(51:340,:));

%%
YYlim=[-0.001 0.015];
figure(123);
subplot(3,8,9)
plotshaded([1:size(Hab1_Posi,1)],[Hab1_Posi_ave+Hab1_Posi_sem Hab1_Posi_ave Hab1_Posi_ave-Hab1_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab1_Nega,1)],[Hab1_Nega_ave+Hab1_Nega_sem Hab1_Nega_ave Hab1_Nega_ave-Hab1_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab1')

subplot(3,8,10)
plotshaded([1:size(Hab2_Posi,1)],[Hab2_Posi_ave+Hab2_Posi_sem Hab2_Posi_ave Hab2_Posi_ave-Hab2_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab2_Nega,1)],[Hab2_Nega_ave+Hab2_Nega_sem Hab2_Nega_ave Hab2_Nega_ave-Hab2_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab2')

subplot(3,8,11)
plotshaded([1:size(FC_Posi,1)],[FC_Posi_ave+FC_Posi_sem FC_Posi_ave FC_Posi_ave-FC_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(FC_Nega,1)],[FC_Nega_ave+FC_Nega_sem FC_Nega_ave FC_Nega_ave-FC_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('FC')

subplot(3,8,12)
plotshaded([1:size(Ext1e_Posi,1)],[Ext1e_Posi_ave+Ext1e_Posi_sem Ext1e_Posi_ave Ext1e_Posi_ave-Ext1e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1e_Nega,1)],[Ext1e_Nega_ave+Ext1e_Nega_sem Ext1e_Nega_ave Ext1e_Nega_ave-Ext1e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial1-5')

subplot(3,8,13)
plotshaded([1:size(Ext1l_Posi,1)],[Ext1l_Posi_ave+Ext1l_Posi_sem Ext1l_Posi_ave Ext1l_Posi_ave-Ext1l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1l_Nega,1)],[Ext1l_Nega_ave+Ext1l_Nega_sem Ext1l_Nega_ave Ext1l_Nega_ave-Ext1l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial21-25')

subplot(3,8,14)
plotshaded([1:size(Ext2e_Posi,1)],[Ext2e_Posi_ave+Ext2e_Posi_sem Ext2e_Posi_ave Ext2e_Posi_ave-Ext2e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2e_Nega,1)],[Ext2e_Nega_ave+Ext2e_Nega_sem Ext2e_Nega_ave Ext2e_Nega_ave-Ext2e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial1-5')

subplot(3,8,15)
plotshaded([1:size(Ext2l_Posi,1)],[Ext2l_Posi_ave+Ext2l_Posi_sem Ext2l_Posi_ave Ext2l_Posi_ave-Ext2l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2l_Nega,1)],[Ext2l_Nega_ave+Ext2l_Nega_sem Ext2l_Nega_ave Ext2l_Nega_ave-Ext2l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial21-25')

subplot(3,8,16)
plotshaded([1:size(Test_Posi,1)],[Test_Posi_ave+Ext2l_Posi_sem Test_Posi_ave Test_Posi_ave-Test_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Test_Nega,1)],[Test_Nega_ave+Test_Nega_sem Test_Nega_ave Test_Nega_ave-Test_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('ExtTest')

%AUC (230624)
figure(124);
subplot(3,9,10)
scatter(0.9,Hab1_Posi_AUC,'k'),hold on,scatter(2.1,Hab1_Nega_AUC,'k')
errorbar(1, mean(Hab1_Posi_AUC), std(Hab1_Posi_AUC)/sqrt(length(Hab1_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab1_Nega_AUC), std(Hab1_Nega_AUC)/sqrt(length(Hab1_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab1  p=' num2str(ranksum(Hab1_Posi_AUC, Hab1_Nega_AUC))])

subplot(3,9,11)
scatter(0.9,Hab2_Posi_AUC,'k'),hold on,scatter(2.1,Hab2_Nega_AUC,'k')
errorbar(1, mean(Hab2_Posi_AUC), std(Hab2_Posi_AUC)/sqrt(length(Hab2_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab2_Nega_AUC), std(Hab2_Nega_AUC)/sqrt(length(Hab2_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab2  p=' num2str(ranksum(Hab2_Posi_AUC, Hab2_Nega_AUC))])

subplot(3,9,12)
scatter(0.9,FCCS_Posi_AUC,'k'),hold on,scatter(2.1,FCCS_Nega_AUC,'k')
errorbar(1, mean(FCCS_Posi_AUC), std(FCCS_Posi_AUC)/sqrt(length(FCCS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCCS_Nega_AUC), std(FCCS_Nega_AUC)/sqrt(length(FCCS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC-CS  p=' num2str(ranksum(FCCS_Posi_AUC, FCCS_Nega_AUC))])

subplot(3,9,13)
scatter(0.9,FCUS_Posi_AUC,'k'),hold on,scatter(2.1,FCUS_Nega_AUC,'k')
errorbar(1, mean(FCUS_Posi_AUC), std(FCUS_Posi_AUC)/sqrt(length(FCUS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCUS_Nega_AUC), std(FCUS_Nega_AUC)/sqrt(length(FCUS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC-US  p=' num2str(ranksum(FCUS_Posi_AUC, FCUS_Nega_AUC))])

subplot(3,9,14)
scatter(0.9,Ext1e_Posi_AUC,'k'),hold on,scatter(2.1,Ext1e_Nega_AUC,'k')
errorbar(1, mean(Ext1e_Posi_AUC), std(Ext1e_Posi_AUC)/sqrt(length(Ext1e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1e_Nega_AUC), std(Ext1e_Nega_AUC)/sqrt(length(Ext1e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1 1-5  p=' num2str(ranksum(Ext1e_Posi_AUC, Ext1e_Nega_AUC))])

subplot(3,9,15)
scatter(0.9,Ext1l_Posi_AUC,'k'),hold on,scatter(2.1,Ext1l_Nega_AUC,'k')
errorbar(1, mean(Ext1l_Posi_AUC), std(Ext1l_Posi_AUC)/sqrt(length(Ext1l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1l_Nega_AUC), std(Ext1l_Nega_AUC)/sqrt(length(Ext1l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1 21-25  p=' num2str(ranksum(Ext1l_Posi_AUC, Ext1l_Nega_AUC))])

subplot(3,9,16)
scatter(0.9,Ext2e_Posi_AUC,'k'),hold on,scatter(2.1,Ext2e_Nega_AUC,'k')
errorbar(1, mean(Ext2e_Posi_AUC), std(Ext2e_Posi_AUC)/sqrt(length(Ext2e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2e_Nega_AUC), std(Ext2e_Nega_AUC)/sqrt(length(Ext2e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2 1-5  p=' num2str(ranksum(Ext2e_Posi_AUC, Ext2e_Nega_AUC))])

subplot(3,9,17)
scatter(0.9,Ext2l_Posi_AUC,'k'),hold on,scatter(2.1,Ext2l_Nega_AUC,'k')
errorbar(1, mean(Ext2l_Posi_AUC), std(Ext2l_Posi_AUC)/sqrt(length(Ext2l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2l_Nega_AUC), std(Ext2l_Nega_AUC)/sqrt(length(Ext2l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2 21-25  p=' num2str(ranksum(Ext2l_Posi_AUC, Ext2l_Nega_AUC))])

subplot(3,9,18)
scatter(0.9,Test_Posi_AUC,'k'),hold on,scatter(2.1,Test_Nega_AUC,'k')
errorbar(1, mean(Test_Posi_AUC), std(Test_Posi_AUC)/sqrt(length(Test_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Test_Nega_AUC), std(Test_Nega_AUC)/sqrt(length(Test_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Test  p=' num2str(ranksum(Test_Posi_AUC, Test_Nega_AUC))])



%% 230625 ave traces (CS responsive during Hab2)

Hab1_Posi=mean(RespMat_Hab1_Posi,3)-RespMatBase_Hab1_Posi;
Hab1_Posi=Hab1_Posi(:,pEexc);
Hab1_Posi_ave=mean(Hab1_Posi,2);
Hab1_Posi_sem=std(Hab1_Posi,[],2)/sqrt(size(Hab1_Posi,2));
Hab1_Nega=mean(RespMat_Hab1_Nega,3)-RespMatBase_Hab1_Nega;
Hab1_Nega=Hab1_Nega(:,nEexc);
Hab1_Nega_ave=mean(Hab1_Nega,2);
Hab1_Nega_sem=std(Hab1_Nega,[],2)/sqrt(size(Hab1_Nega,2));

Hab2_Posi=mean(RespMat_Hab2_Posi,3)-RespMatBase_Hab2_Posi;
Hab2_Posi=Hab2_Posi(:,pEexc);
Hab2_Posi_ave=mean(Hab2_Posi,2);
Hab2_Posi_sem=std(Hab2_Posi,[],2)/sqrt(size(Hab2_Posi,2));
Hab2_Nega=mean(RespMat_Hab2_Nega,3)-RespMatBase_Hab2_Nega;
Hab2_Nega=Hab2_Nega(:,nEexc);
Hab2_Nega_ave=mean(Hab2_Nega,2);
Hab2_Nega_sem=std(Hab2_Nega,[],2)/sqrt(size(Hab2_Nega,2));

FC_Posi=mean(RespMat_FC_Posi,3)-RespMatBase_FC_Posi;
FC_Posi=FC_Posi(:,pEexc);
FC_Posi_ave=mean(FC_Posi,2);
FC_Posi_sem=std(FC_Posi,[],2)/sqrt(size(FC_Posi,2));
FC_Nega=mean(RespMat_FC_Nega,3)-RespMatBase_FC_Nega;
FC_Nega=FC_Nega(:,nEexc);
FC_Nega_ave=mean(FC_Nega,2);
FC_Nega_sem=std(FC_Nega,[],2)/sqrt(size(FC_Nega,2));

Ext1e_Posi=mean(RespMat_Ext1_Posi(:,:,1:5),3)-RespMatBase_Ext1e_Posi;
Ext1e_Posi=Ext1e_Posi(:,pEexc);
Ext1e_Posi_ave=mean(Ext1e_Posi,2);
Ext1e_Posi_sem=std(Ext1e_Posi,[],2)/sqrt(size(Ext1e_Posi,2));
Ext1e_Nega=mean(RespMat_Ext1_Nega(:,:,1:5),3)-RespMatBase_Ext1e_Nega;
Ext1e_Nega=Ext1e_Nega(:,nEexc);
Ext1e_Nega_ave=mean(Ext1e_Nega,2);
Ext1e_Nega_sem=std(Ext1e_Nega,[],2)/sqrt(size(Ext1e_Nega,2));

Ext1l_Posi=mean(RespMat_Ext1_Posi(:,:,21:25),3)-RespMatBase_Ext1l_Posi;
Ext1l_Posi=Ext1l_Posi(:,pEexc);
Ext1l_Posi_ave=mean(Ext1l_Posi,2);
Ext1l_Posi_sem=std(Ext1l_Posi,[],2)/sqrt(size(Ext1l_Posi,2));
Ext1l_Nega=mean(RespMat_Ext1_Nega(:,:,21:25),3)-RespMatBase_Ext1l_Nega;
Ext1l_Nega=Ext1l_Nega(:,nEexc);
Ext1l_Nega_ave=mean(Ext1l_Nega,2);
Ext1l_Nega_sem=std(Ext1l_Nega,[],2)/sqrt(size(Ext1l_Nega,2));

Ext2e_Posi=mean(RespMat_Ext2_Posi(:,:,1:5),3)-RespMatBase_Ext2e_Posi;
Ext2e_Posi=Ext2e_Posi(:,pEexc);
Ext2e_Posi_ave=mean(Ext2e_Posi,2);
Ext2e_Posi_sem=std(Ext2e_Posi,[],2)/sqrt(size(Ext2e_Posi,2));
Ext2e_Nega=mean(RespMat_Ext2_Nega(:,:,1:5),3)-RespMatBase_Ext2e_Nega;
Ext2e_Nega=Ext2e_Nega(:,nEexc);
Ext2e_Nega_ave=mean(Ext2e_Nega,2);
Ext2e_Nega_sem=std(Ext2e_Nega,[],2)/sqrt(size(Ext2e_Nega,2));

Ext2l_Posi=mean(RespMat_Ext2_Posi(:,:,21:25),3)-RespMatBase_Ext2l_Posi;
Ext2l_Posi=Ext2l_Posi(:,pEexc);
Ext2l_Posi_ave=mean(Ext2l_Posi,2);
Ext2l_Posi_sem=std(Ext2l_Posi,[],2)/sqrt(size(Ext2l_Posi,2));
Ext2l_Nega=mean(RespMat_Ext2_Nega(:,:,21:25),3)-RespMatBase_Ext2l_Nega;
Ext2l_Nega=Ext2l_Nega(:,nEexc);
Ext2l_Nega_ave=mean(Ext2l_Nega,2);
Ext2l_Nega_sem=std(Ext2l_Nega,[],2)/sqrt(size(Ext2l_Nega,2));

Test_Posi=mean(RespMat_Test_Posi,3)-RespMatBase_Test_Posi;
Test_Posi=Test_Posi(:,pEexc);
Test_Posi_ave=mean(Test_Posi,2);
Test_Posi_sem=std(Test_Posi,[],2)/sqrt(size(Test_Posi,2));
Test_Nega=mean(RespMat_Test_Nega,3)-RespMatBase_Test_Nega;
Test_Nega=Test_Nega(:,nEexc);
Test_Nega_ave=mean(Test_Nega,2);
Test_Nega_sem=std(Test_Nega,[],2)/sqrt(size(Test_Nega,2));

%AUC (230624)
Hab1_Posi_AUC=mean(Hab1_Posi(51:340,:));
Hab2_Posi_AUC=mean(Hab2_Posi(51:340,:));
FCCS_Posi_AUC=mean(FC_Posi(51:340,:));
FCUS_Posi_AUC=mean(FC_Posi(340:380,:));
Ext1e_Posi_AUC=mean(Ext1e_Posi(51:340,:));
Ext1l_Posi_AUC=mean(Ext1l_Posi(51:340,:));
Ext2e_Posi_AUC=mean(Ext2e_Posi(51:340,:));
Ext2l_Posi_AUC=mean(Ext2l_Posi(51:340,:));
Test_Posi_AUC=mean(Test_Posi(51:340,:));

Hab1_Nega_AUC=mean(Hab1_Nega(51:340,:));
Hab2_Nega_AUC=mean(Hab2_Nega(51:340,:));
FCCS_Nega_AUC=mean(FC_Nega(51:340,:));
FCUS_Nega_AUC=mean(FC_Nega(340:380,:));
Ext1e_Nega_AUC=mean(Ext1e_Nega(51:340,:));
Ext1l_Nega_AUC=mean(Ext1l_Nega(51:340,:));
Ext2e_Nega_AUC=mean(Ext2e_Nega(51:340,:));
Ext2l_Nega_AUC=mean(Ext2l_Nega(51:340,:));
Test_Nega_AUC=mean(Test_Nega(51:340,:));

%%
YYlim=[-0.001 0.015];

figure(123);
subplot(3,8,17)
plotshaded([1:size(Hab1_Posi,1)],[Hab1_Posi_ave+Hab1_Posi_sem Hab1_Posi_ave Hab1_Posi_ave-Hab1_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab1_Nega,1)],[Hab1_Nega_ave+Hab1_Nega_sem Hab1_Nega_ave Hab1_Nega_ave-Hab1_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab1')

subplot(3,8,18)
plotshaded([1:size(Hab2_Posi,1)],[Hab2_Posi_ave+Hab2_Posi_sem Hab2_Posi_ave Hab2_Posi_ave-Hab2_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Hab2_Nega,1)],[Hab2_Nega_ave+Hab2_Nega_sem Hab2_Nega_ave Hab2_Nega_ave-Hab2_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Hab2')

subplot(3,8,19)
plotshaded([1:size(FC_Posi,1)],[FC_Posi_ave+FC_Posi_sem FC_Posi_ave FC_Posi_ave-FC_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(FC_Nega,1)],[FC_Nega_ave+FC_Nega_sem FC_Nega_ave FC_Nega_ave-FC_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('FC')

subplot(3,8,20)
plotshaded([1:size(Ext1e_Posi,1)],[Ext1e_Posi_ave+Ext1e_Posi_sem Ext1e_Posi_ave Ext1e_Posi_ave-Ext1e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1e_Nega,1)],[Ext1e_Nega_ave+Ext1e_Nega_sem Ext1e_Nega_ave Ext1e_Nega_ave-Ext1e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial1 5')

subplot(3,8,21)
plotshaded([1:size(Ext1l_Posi,1)],[Ext1l_Posi_ave+Ext1l_Posi_sem Ext1l_Posi_ave Ext1l_Posi_ave-Ext1l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext1l_Nega,1)],[Ext1l_Nega_ave+Ext1l_Nega_sem Ext1l_Nega_ave Ext1l_Nega_ave-Ext1l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext1: trial21 25')

subplot(3,8,22)
plotshaded([1:size(Ext2e_Posi,1)],[Ext2e_Posi_ave+Ext2e_Posi_sem Ext2e_Posi_ave Ext2e_Posi_ave-Ext2e_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2e_Nega,1)],[Ext2e_Nega_ave+Ext2e_Nega_sem Ext2e_Nega_ave Ext2e_Nega_ave-Ext2e_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial1 5')

subplot(3,8,23)
plotshaded([1:size(Ext2l_Posi,1)],[Ext2l_Posi_ave+Ext2l_Posi_sem Ext2l_Posi_ave Ext2l_Posi_ave-Ext2l_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Ext2l_Nega,1)],[Ext2l_Nega_ave+Ext2l_Nega_sem Ext2l_Nega_ave Ext2l_Nega_ave-Ext2l_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('Ext2: trial21 25')

subplot(3,8,24)
plotshaded([1:size(Test_Posi,1)],[Test_Posi_ave+Ext2l_Posi_sem Test_Posi_ave Test_Posi_ave-Test_Posi_sem],[0.9 0 0],1),hold on
plotshaded([1:size(Test_Nega,1)],[Test_Nega_ave+Test_Nega_sem Test_Nega_ave Test_Nega_ave-Test_Nega_sem],[0 0 0.9],1)
plot([0 size(Hab1_Posi,1)],[0 0],'k--')
plot([50 50],[YYlim],'k--')
plot([340 340],[YYlim],'k--')
ylim([YYlim])
title('ExtTest')


%AUC (230624)
figure(124);
subplot(3,9,19)
scatter(0.9,Hab1_Posi_AUC,'k'),hold on,scatter(2.1,Hab1_Nega_AUC,'k')
errorbar(1, mean(Hab1_Posi_AUC), std(Hab1_Posi_AUC)/sqrt(length(Hab1_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab1_Nega_AUC), std(Hab1_Nega_AUC)/sqrt(length(Hab1_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab1  p=' num2str(ranksum(Hab1_Posi_AUC, Hab1_Nega_AUC))])

subplot(3,9,20)
scatter(0.9,Hab2_Posi_AUC,'k'),hold on,scatter(2.1,Hab2_Nega_AUC,'k')
errorbar(1, mean(Hab2_Posi_AUC), std(Hab2_Posi_AUC)/sqrt(length(Hab2_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Hab2_Nega_AUC), std(Hab2_Nega_AUC)/sqrt(length(Hab2_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Hab2  p=' num2str(ranksum(Hab2_Posi_AUC, Hab2_Nega_AUC))])

subplot(3,9,21)
scatter(0.9,FCCS_Posi_AUC,'k'),hold on,scatter(2.1,FCCS_Nega_AUC,'k')
errorbar(1, mean(FCCS_Posi_AUC), std(FCCS_Posi_AUC)/sqrt(length(FCCS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCCS_Nega_AUC), std(FCCS_Nega_AUC)/sqrt(length(FCCS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC CS  p=' num2str(ranksum(FCCS_Posi_AUC, FCCS_Nega_AUC))])

subplot(3,9,22)
scatter(0.9,FCUS_Posi_AUC,'k'),hold on,scatter(2.1,FCUS_Nega_AUC,'k')
errorbar(1, mean(FCUS_Posi_AUC), std(FCUS_Posi_AUC)/sqrt(length(FCUS_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(FCUS_Nega_AUC), std(FCUS_Nega_AUC)/sqrt(length(FCUS_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['FC US  p=' num2str(ranksum(FCUS_Posi_AUC, FCUS_Nega_AUC))])

subplot(3,9,23)
scatter(0.9,Ext1e_Posi_AUC,'k'),hold on,scatter(2.1,Ext1e_Nega_AUC,'k')
errorbar(1, mean(Ext1e_Posi_AUC), std(Ext1e_Posi_AUC)/sqrt(length(Ext1e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1e_Nega_AUC), std(Ext1e_Nega_AUC)/sqrt(length(Ext1e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1 1-5  p=' num2str(ranksum(Ext1e_Posi_AUC, Ext1e_Nega_AUC))])

subplot(3,9,24)
scatter(0.9,Ext1l_Posi_AUC,'k'),hold on,scatter(2.1,Ext1l_Nega_AUC,'k')
errorbar(1, mean(Ext1l_Posi_AUC), std(Ext1l_Posi_AUC)/sqrt(length(Ext1l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext1l_Nega_AUC), std(Ext1l_Nega_AUC)/sqrt(length(Ext1l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext1 1-5  p=' num2str(ranksum(Ext1l_Posi_AUC, Ext1l_Nega_AUC))])

subplot(3,9,25)
scatter(0.9,Ext2e_Posi_AUC,'k'),hold on,scatter(2.1,Ext2e_Nega_AUC,'k')
errorbar(1, mean(Ext2e_Posi_AUC), std(Ext2e_Posi_AUC)/sqrt(length(Ext2e_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2e_Nega_AUC), std(Ext2e_Nega_AUC)/sqrt(length(Ext2e_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2 1-5  p=' num2str(ranksum(Ext2e_Posi_AUC, Ext2e_Nega_AUC))])

subplot(3,9,26)
scatter(0.9,Ext2l_Posi_AUC,'k'),hold on,scatter(2.1,Ext2l_Nega_AUC,'k')
errorbar(1, mean(Ext2l_Posi_AUC), std(Ext2l_Posi_AUC)/sqrt(length(Ext2l_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Ext2l_Nega_AUC), std(Ext2l_Nega_AUC)/sqrt(length(Ext2l_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Ext2 21-25  p=' num2str(ranksum(Ext2l_Posi_AUC, Ext2l_Nega_AUC))])

subplot(3,9,27)
scatter(0.9,Test_Posi_AUC,'k'),hold on,scatter(2.1,Test_Nega_AUC,'k')
errorbar(1, mean(Test_Posi_AUC), std(Test_Posi_AUC)/sqrt(length(Test_Posi_AUC)),'r',LineWidth=2)
errorbar(2, mean(Test_Nega_AUC), std(Test_Nega_AUC)/sqrt(length(Test_Nega_AUC)),'b',LineWidth=2)
plot([0.5 2.5],[0 0],'--'),xlim([0.5 2.5]),ylim([-0.005 0.02])
title(['Test  p=' num2str(ranksum(Test_Posi_AUC, Test_Nega_AUC))])