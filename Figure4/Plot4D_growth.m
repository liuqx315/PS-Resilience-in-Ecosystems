clear all;
close all;
cgray = gray(10);

SDF=load('Data/SDF_data.mat');
DDA=load('Data/DDA_data.mat');
%% ==========================
dat_SDF=SDF.A;
dat_DDA=DDA.A;
img1=dat_SDF(:,:,50); % change here for SDF or DDA analysis
img2=dat_SDF(:,:,100);
Rmin = 2;
Rmax = 25;

BWimg1=img1>mean(img1(:)); % imbinarize(img1);
BWimg2=img2>mean(img2(:));

stats1=regionprops('table',BWimg1,'Centroid','Area','Circularity','EquivDiameter',...
    'Orientation','Eccentricity','MajorAxisLength','MinorAxisLength'); % 'all'

index1=1:size(stats1,1);
centroids1 = [index1' stats1.Area cat(1,stats1.Centroid) stats1.Circularity stats1.EquivDiameter stats1.Orientation];

stats2=regionprops('table',BWimg2,'Centroid','Area','Circularity','EquivDiameter',...
    'Orientation','Eccentricity','MajorAxisLength','MinorAxisLength'); % 'all'
index2=1:size(stats2,1);
centroids2 = [index2' stats2.Area cat(1,stats2.Centroid) stats2.Circularity stats2.EquivDiameter stats2.Orientation];


for k=1:size(centroids1,1)
    
    distMin=100000;
    index=1;
    for j=1:size(centroids2,1)
     dist=sqrt((centroids1(k,3)-centroids2(j,3))^2+(centroids1(k,4)-centroids2(j,4))^2);
     if distMin>dist
         distMin=dist;
         indx=j;
     end
    end
     dat(k,:)=[centroids1(k,:) centroids2(indx,:)];
end


xlswrite("SDF_50",centroids1)
xlswrite("SDF_100",centroids2)
xlswrite("SDF_50vs100",dat)

Growth1=dat(:,9)-dat(:,2);
Growth2=dat(:,13)-dat(:,6);

figure('position', [100,10,2000,2000])
plot(dat(:,2),Growth1,'o');
% hold on
% text(dat(:,2)-8,Growth,num2str(dat(:,1)),'Color','red','FontSize',10)
% text(dat(:,2)-8,Growth-10,num2str(dat(:,8)),'Color',cgray(7,:),'FontSize',10)

% figure('position', [100,10,2000,2000])
% t = tiledlayout(1,1,'TileSpacing','Compact');
% t.Padding = 'compact';

imshow(BWimg1)
hold on
text(dat(:,3)-8,dat(:,4),num2str(dat(:,1)),'Color','red','FontSize',5)
saveas(gcf,'SDF_50.png')
