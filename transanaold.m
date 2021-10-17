clear all 
clc

% dynamic 
expdate = 20170308%%%%%% experiment date
phstart =1;  %%%%%%%%%%%%%
totalslice =20; %%%%
slice = totalslice; %%%%%%%%%  % display how many slices or measurements; 
dispx=5; dispy=9;    %%%%       % display how main one figure;
interval=1;          %%%%
brightness = 0.01;
quater = 4; 
HIFU_start=5;
HIFU_end=15;
Venc=6; %% % cm/s
chamberdiameter = 2.2; %%%%% mm
chamberradius = chamberdiameter/2; %%%%%%% mm
conc = 0.1; %%%%%% concerntration in %
magnum = 132
phsnum = 133
xlsfile = ['ROI_', num2str(conc),'%','.xls'];
histfile = ['hist', num2str(expdate),'_',num2str(conc),'%','.xls'];
dicomfile=[num2str(expdate),'_HLL_RK_',num2str(expdate),'_s0',num2str(phsnum),'_i00%03d'];
dicomfile1=[num2str(expdate),'_HLL_RK_',num2str(expdate),'_s0',num2str(magnum),'_i00%03d'];
% ROI_image_focus_mag = dicomread('20170308_HLL_RK_20170308_s0144_i00003.dcm');
% ROI_image_focus_phs = dicomread('20170308_HLL_RK_20170308_s0145_i00013.dcm');
% ROI_image_bandwith = dicomread('20170308_HLL_RK_20170308_s0145_i00013.dcm');
%% ------ read file ------%
% dynamic images;
count=0;
for ii= phstart:interval:phstart+slice-1
    count=count+1;
    Img(:,:,count) = dicomread(strcat(sprintf(dicomfile,ii),'.dcm'));
    if count < 2
        temp_info=dicominfo(strcat(sprintf(dicomfile,ii),'.dcm'));
        SliceThickness=temp_info.SliceThickness;
        %SliceGap=temp_info.SpacingBetweenSlices-SliceThickness;
        TR = temp_info.RepetitionTime;
        TE = temp_info.EchoTime;
        EchoTranLen = temp_info.EchoTrainLength;
        PixelSize=temp_info.PixelSpacing;
        matrix_x = temp_info.Rows;
        matrix_y = temp_info.Columns;
        NEX=temp_info.NumberOfAverages;
    end
end
count=0;
for ii= phstart:interval:phstart+slice-1
    count=count+1;
    Img_mag(:,:,count) = dicomread(strcat(sprintf(dicomfile1,ii),'.dcm'));
end
FOV_x = matrix_y*PixelSize(1);
FOV_y = matrix_x*PixelSize(2);
time = double((TR*matrix_y*NEX))*([1:size(Img,3)]-1)./1000;  % ms=>s;
TIME = reshape(time,[],1);
chamber_in_pixel = chamberradius/PixelSize(1);
Img=double(Img);
 Img_mag=double( Img_mag);
baseline = mean(Img(:,:,2:4),3);
baseline2 = repmat(baseline,[1,1,20]);
nor_Img = ((Img-baseline2)./baseline2)*100;
nor_Img = -nor_Img;
Img = (Img-2048)/2048*Venc;
Img = -Img;
%% ------ selection of stationary tissue ------%
presl = 3;
postsl = 17;
intnum = postsl-presl+1;

Img_mean = mean(Img(:,:,presl:postsl),3);%把所有images平均
temp5 = repmat(Img_mean,[1,1,intnum]);%
temp0 = Img(:,:,presl:postsl)-temp5;

for t = 1:intnum
    for jj = 1:matrix_x
        for ii = 1:matrix_y
        Img_square(jj,ii,t) = temp0(jj,ii,t)^2;
        end
    end
end

Img_sum = sum(Img_square(:,:,1:intnum),3);


 for jj = 1:matrix_x
     for ii = 1:matrix_y
        Img_std(jj,ii) = ((Img_sum(jj,ii)/(intnum)))^(1/2);
     end
 end

Velocity_map = Img(:,:,14);

figure(11)
subplot(1,2,1)
imshow(Velocity_map,[-2 2]);colormap jet
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,2,2)

imshow(Img_std,[0 0.2]);colormap jet
set(gcf,'outerposition',get(0,'screensize'));
%% 
k = max(Img_mag(:,:,3));
kk = max(k(:));
colormin = 0; colormax = kk;

rr = 55; ll = 75; uu = 55; dd = 75;
figure(111)

imshow(Img_mag(:,:,3),[colormin colormax]);colormap gray
title('please select ROI center')
set(gcf,'outerposition',get(0,'screensize'));axis([rr ll uu dd]);
hold on
R0=ginput(1);
R0 = round(R0);
hold off
theta = linspace(0,2*pi,100);
for ii = 1:100
    X(ii) = R0(1)+chamber_in_pixel*sin(theta(ii));
    Y(ii) = R0(2)+chamber_in_pixel*cos(theta(ii));
end
X = reshape(X,[],1);
Y = reshape(Y,[],1);
roi(:,:,1) = roipoly(Img_std,X,Y);
temp = [find(roi)];
SI_temporal_STD = Img_std([temp]);
save R0 'R0'
%% different phase% 
figure(600)
imshow(Img(:,:,presl),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

figure(700)
imshow(Img(:,:,12),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

figure(800)
imshow(Img(:,:,postsl),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

saveas(figure(600),'MB_pre.png');
saveas(figure(700),'MB_FUS.png');
saveas(figure(800),'MB_post.png');
%% different phase mag% 
figure(601)
imshow( Img_mag(:,:,presl),[colormin colormax]);colormap gray;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

figure(701)
imshow( Img_mag(:,:,12),[colormin colormax]);colormap gray;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

figure(801)
imshow( Img_mag(:,:,postsl),[colormin colormax]);colormap gray;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);

saveas(figure(600),'MB_magpre.png');
saveas(figure(700),'MB_magFUS.png');
saveas(figure(800),'MB_magpost.png');
%%  histogram %%
temp_STD_mean = mean(SI_temporal_STD);
temp_STD_std = std(SI_temporal_STD);
temp_STD_sum = sum(SI_temporal_STD);
temp_STD_median = median(SI_temporal_STD);
temp_STD_skew = skewness(SI_temporal_STD);
temp_STD_kurt = kurtosis(SI_temporal_STD);
temp_STD_prc = [prctile(SI_temporal_STD,10);prctile(SI_temporal_STD,25);prctile(SI_temporal_STD,50);prctile(SI_temporal_STD,75);prctile(SI_temporal_STD,90)];
index_STD = [ temp_STD_mean ; temp_STD_std ;temp_STD_sum;temp_STD_median ;temp_STD_skew ;temp_STD_kurt ;temp_STD_prc];
index_name = {'MEAN' ;'STD' ;'SUM' ;'MEDIAN' ;'SKEWNESS' ;'KURTOSIS' ;'PRC 10%' ;'PRC 25%' ;'PRC 50%' ;'PRC 75%' ;'PRC 90%'};
%%
figure(999)
nbins=(max(SI_temporal_STD)-min(SI_temporal_STD))/0.01
hist(SI_temporal_STD,nbins);
xlim([0 0.5]);ylim([0 40]);
h = findobj(gca,'type','patch');
set(gca,'FontWeight','bold','fontsize',16,'linewidth',3)
set(h,'FaceColor','r','EdgeColor','k','linewidth',3)
xlabel('STD','fontsize',16);ylabel('# of pixel','fontsize',16);
saveas(figure(999),'hist_5-25','png');
%% plot on image %%
k = max(Img_mag(:,:,3));
kk = max(k(:));
colormin = 0; colormax = kk;
figure(55555)
imshow(Img_mag(:,:,3),[colormin colormax]);colormap gray
set(gcf,'outerposition',[0 40 1275 717]);axis([rr ll uu dd]);
hold on
plot(X,Y,'y-','linewidth',3)
hold off


figure(555555)
imshow(Img_std,[0 0.4]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
plot(X,Y,'y-','linewidth',3)
hold off

saveas(figure(555555),['STD_phs_',num2str(presl),'-',num2str(postsl),'.png']);
saveas(figure(55555),['STD_mag_',num2str(presl),'-',num2str(postsl),'.png']);
%% find MB %%
ddd = roi(:,:,1).*Img_std;
[row1 col1]=find(ddd);

MB=prctile(SI_temporal_STD,90);
[row col]=find(ddd>MB);
figure(5555555)
imshow(ddd,[0 0.4]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
plot(col,row,'g*','linewidth',5)
hold off
saveas(figure(5555555),'MB-90.png');
%%
MBselect = [row,col];
[idxMB,ctrMB] = kmeans(MBselect,2);
velocity = [MBselect(idxMB==2,1),MBselect(idxMB==2,2)];
mbposition = [MBselect(idxMB==1,1),MBselect(idxMB==1,2)];

V = zeros(128,128);
for jj = 1:size(velocity,1)
  V(velocity(jj,1),velocity(jj,2))=1;
end
for ii = presl:postsl
   aaa = Img(:,:,ii);
   velo_part = find(aaa.*V~=0);
   mean_velo(ii) = mean(aaa(velo_part));
end

B = zeros(128,128);
for jj = 1:size(mbposition,1)
  B(mbposition(jj,1),mbposition(jj,2))=1;
end
for ii = presl:postsl
   aaa = Img(:,:,ii);
   mb_part = find(aaa.*B~=0);
   mean_mb(ii) = mean(aaa(mb_part));
end

mean_velo = reshape(mean_velo,[],1);
mean_mb = reshape(mean_mb,[],1);

figure(6666)
plot(mean_velo,'b','linewidth',3)
title('mean _ velocity','FontWeight','bold','fontsize',10)
% xlim([5 25]);ylim([0 1]);
set(gca,'FontWeight','bold','fontsize',16,'linewidth',3)

figure(66666)
plot(mean_mb,'r','linewidth',3)
title('mean _ mb','FontWeight','bold','fontsize',10)
% xlim([5 25]);ylim([0 0.5]);
set(gca,'FontWeight','bold','fontsize',16,'linewidth',3)
%% normalize
for ii = presl:postsl
   ccc = nor_Img(:,:,ii);
   nor_velo_part = find(ccc.*V~=0);
   nor_mean_velo(ii) = mean(ccc(velo_part));
end
nor_mean_velo = reshape(nor_mean_velo,[],1);
for ii = presl:postsl
   ccc = nor_Img(:,:,ii);
   nor_mb_part = find(ccc.*B~=0);
   nor_mean_mb(ii) = mean(ccc(mb_part));
end
nor_mean_mb = reshape(nor_mean_mb,[],1);
figure(7777)
plot(nor_mean_velo ,'b','linewidth',3)
title('mean _ velocity','FontWeight','bold','fontsize',10)
% xlim([5 25]);ylim([0 1]);
set(gca,'FontWeight','bold','fontsize',16,'linewidth',3)
figure(77777)
plot(nor_mean_mb ,'r','linewidth',3)
title('mean _ velocity','FontWeight','bold','fontsize',10)
% xlim([5 25]);ylim([0 1]);
set(gca,'FontWeight','bold','fontsize',16,'linewidth',3)
%% % kmeans
figure(500)
imshow(ddd,[0 0.4]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
plot(MBselect(idxMB==1,2),MBselect(idxMB==1,1),'k*','linewidth',5)
plot(MBselect(idxMB==2,2),MBselect(idxMB==2,1),'r*','linewidth',5)
hold off
saveas(figure(500),'MBkmeans_90.png');
%% kmeans on velomap % 
figure(602)
imshow(Img(:,:,presl),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
% plot(col,row,'b*','linewidth',5)
plot(MBselect(idxMB==1,2),MBselect(idxMB==1,1),'k*','linewidth',5)
plot(MBselect(idxMB==2,2),MBselect(idxMB==2,1),'r*','linewidth',5)
hold off
figure(702)
imshow(Img(:,:,12),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
% plot(col,row,'b*','linewidth',5)
plot(MBselect(idxMB==1,2),MBselect(idxMB==1,1),'k*','linewidth',5)
plot(MBselect(idxMB==2,2),MBselect(idxMB==2,1),'r*','linewidth',5)
hold off
figure(802)
imshow(Img(:,:,postsl),[-1 1.5]);colormap jet;colorbar 
set(gcf,'outerposition',[0 40 1275 717]); axis([rr ll uu dd]);
hold on
% plot(col,row,'b*','linewidth',5)
plot(MBselect(idxMB==1,2),MBselect(idxMB==1,1),'k*','linewidth',5)
plot(MBselect(idxMB==2,2),MBselect(idxMB==2,1),'r*','linewidth',5)
hold off
saveas(figure(600),'MBkmeans_pre.png');
saveas(figure(700),'MBkmeans_FUS.png');
saveas(figure(800),'MBkmeans_post.png');
%% write to excel %%
xlswrite(xlsfile,R0,'chamber center','B2');
xlswrite(xlsfile,SI_temporal_STD,'SI_temporal_STD','C4');
xlswrite(xlsfile,index_name,'statistic state','A2');
xlswrite(xlsfile,index_STD,'statistic state','B2');
xlswrite(xlsfile,TIME,'velocity','A2');
xlswrite(xlsfile,mean_velo,'velocity','B2');
xlswrite(xlsfile,mean_mb,'velocity','C2');
xlswrite(xlsfile,TIME,'nor_velocity','A2');
xlswrite(xlsfile,nor_mean_velo,'nor_velocity','B2');
xlswrite(xlsfile,nor_mean_mb,'nor_velocity','C2');
% 
% close all