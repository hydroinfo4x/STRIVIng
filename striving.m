%%INFO
%Plotting WorldClim Data processed by Vitali Diaz Mercado
%-for the baseline and two periods: 2050 (2041-2060) and 2070 (2061-2080)
%-for 4 RCPs
%-plotting 20 maps: baseline + 19 GCM with climate scanarios
%-plotting monthly hietograms and temperatures to research
% Changes in interannual variability (see poster)
% 
%Check the updates at https://github.com/hydroinfo4x/STRIVIng

%%
clear all;
close all;
clc;

%%
%%Computing fi
xmin=-180;
xmax=180;
ymin=-60;
ymax=90;
cols=2160;
rows=900;

yres=(ymax-ymin)/rows;
lat_deg=(ymin+yres/2):yres:(ymax-yres/2);
lat_deg=sort(lat_deg,'descend');
lat_mat=zeros(rows,cols);
for jj=1:cols
    lat_mat(:,jj)=lat_deg(:)*pi()/180; %fi
end

%%Computing delta and dr
J=1:365;
D=80;
d=(23.45*pi()/180)*sin(2*pi()*(J-D)/365); %delta
dr=1+0.033*cos(2*pi()*J/365); %dr

days_month(1,:)=[1,32,60,91,121,152,182,213,244,274,305,335];
days_month(2,:)=[31,59,90,120,151,181,212,243,273,304,334,365];
d_vec=[];
for mm=1:12
    d_vec(mm)=mean(d(days_month(1,mm):days_month(2,mm)));
end
dr_vec=[];
for mm=1:12
    dr_vec(mm)=mean(dr(days_month(1,mm):days_month(2,mm)));
end

%%Computing ws
ws_month=[];
for mm=1:12
    ws_month(:,:,mm)=-tan(lat_mat)*tan(d_vec(mm));
end


%%  NAMES
%(1.1)Model Names 
model_label=cell(20,1);
model_label{1,1}='Current (~1950-2000)';
model_label{2,1}='ACCESS1-0';
model_label{3,1}='BCC-CSM1-1';
model_label{4,1}='CCSM4';
model_label{5,1}='CESM1-CAM5-1-FV2';
model_label{6,1}='CNRM-CM5';
model_label{7,1}='GFDL-ESM2G';
model_label{8,1}='GFDL-CM3';
model_label{9,1}='GISS-E2-R';
model_label{10,1}='HadGEM2-AO';
model_label{11,1}='HadGEM2-ES';
model_label{12,1}='HadGEM2-CC';
model_label{13,1}='INMCM4';
model_label{14,1}='IPSL-CM5A-LR';
model_label{15,1}='MIROC5';
model_label{16,1}='MRI-CGCM3';
model_label{17,1}='MIROC-ESM-CHEM';
model_label{18,1}='MPI-ESM-LR';
model_label{19,1}='MIROC-ESM';
model_label{20,1}='NorESM1-M';

%Variable names
variable_label=cell(6,1);
variable_label{1,1}='Precipitation';
variable_label{2,1}='Mean temperature';
variable_label{3,1}='Minimum temperature';
variable_label{4,1}='Maximum temperature';
variable_label{5,1}='Precipitation';
variable_label{6,1}='Temperature';

%Period names
period_label=cell(2,1);
period_label{1,1}='2050 (2041-2060)'; %'2041-2060'
period_label{2,1}='2070 (2061-2080)'; %'2061-2080'

%Path names
path_label=cell(4,1);
path_label{1,1}='RCP 2.6';
path_label{2,1}='RCP 4.5';
path_label{3,1}='RCP 6.0';
path_label{4,1}='RCP 8.5';

%Subvariable
subvariable_label=cell(12,6);
%Precipitation,temperatures
sublabel_units=' (mm)';
for jj=1:4
    if jj>1
        sublabel_units=' (°C)';
    end
    subvariable_label{1,jj}=['Jan',sublabel_units];
    subvariable_label{2,jj}=['Feb',sublabel_units];
    subvariable_label{3,jj}=['Mar',sublabel_units];
    subvariable_label{4,jj}=['Apr',sublabel_units];
    subvariable_label{5,jj}=['May',sublabel_units];
    subvariable_label{6,jj}=['Jun',sublabel_units];
    subvariable_label{7,jj}=['Jul',sublabel_units];
    subvariable_label{8,jj}=['Aug',sublabel_units];
    subvariable_label{9,jj}=['Sep',sublabel_units];
    subvariable_label{10,jj}=['Oct',sublabel_units];
    subvariable_label{11,jj}=['Nov',sublabel_units];
    subvariable_label{12,jj}=['Dec',sublabel_units];
end

%BIOs
%Precipitation
subvariable_label{1,5}='Annual (mm)';
subvariable_label{2,5}='Wettest Month (mm)';
subvariable_label{3,5}='Driest Month (mm)';
subvariable_label{4,5}='Seasonality';
subvariable_label{5,5}='Wettest Quarter (mm)';
subvariable_label{6,5}='Driest Quarter (mm)';
subvariable_label{7,5}='Warmer Quarter (mm)';
subvariable_label{8,5}='Coldest Quarter (mm)';

%Precipitation
subvariable_label{1,6}='Annual mean (°C)';
subvariable_label{2,6}='Maximum of Warmest Month (°C)';
subvariable_label{3,6}='Minimum of Coldest Month (°C)';
subvariable_label{4,6}='Seasonality (*100)';
subvariable_label{5,6}='Mean of Wettest Quarter (°C)';
subvariable_label{6,6}='Mean of Driest Quarter (°C)';
subvariable_label{7,6}='Mean of Warmer Quarter (°C)';
subvariable_label{8,6}='Mean of Coldest Quarter (°C)';
subvariable_label{9,6}='Mean Diurnal Range (°C)';
subvariable_label{10,6}='Annual Range (°C)';
subvariable_label{11,6}='Isothermality';


%%
%***********************
%(1.2)CONFIGURATION: MATs
folderbigMAT='D:\02_DATA_D\WorldClim\CMIP5\MATs\';
suffix='wc';

%Choosing variable
variable_opc = 6; %val: 1)pr, 2)tm, 3)tn, 4)tx, 5)bp, 6)bt
id_variable=1;  %1-12) 

%Only one
period_opc=1; %period: 1)50, 2)70
path_opc=2; %path: 1)26, 2)45, 3)60, 4)85

%Batch
% pi=1; %period: 1)50, 2)70
% pf=2;
% ri=1; %path: 1)26, 2)45, 3)60, 4)85
% rf=4;

% %*****GLOBE*****
% nameMask='wcmask10m';
% xmin=1;  %GLOBALxmin=1; 
% xmax=cols;  %xmax=cols;
% ymin=1;  %ymin=1;
% ymax=rows;  %ymax=rows;
% mask_Area=510072000-14000000;  %Globe: 510072000 km2, Antartica: 14,000,000 km2
% delta_axis=30;  %10 degrees

% %*****REGION*****
% %%Latin-America
% xmin=360;
% xmax=887;
% ymin=310;
% ymax=880;
% mask_Area=510072000;  %Globe: 510,072,000 km2
% delta_axis=20;  %10 degrees

% %%Mexico
% %Mexico
% nameMask='wcmaskMexico10m';
% xmin=360;
% xmax=565;
% ymin=310;
% ymax=460;
% mask_Area=1972550;  %Globe: 510,072,000 km2
% delta_axis=10;  %10 degrees

%%Dominican Republic
nameMask='wcmaskDominicanRepu10m';
xmin=648;
xmax=673;
ymin=419;
ymax=437;
mask_Area=48670;  %Globe: 510,072,000 km2
delta_axis=1;  %10 degrees
% 
% %%Amazon Basin
% nameMask='wcmask_Amazon_Basin';
% xmin=600;
% xmax=781;
% ymin=506;
% ymax=667;
% mask_Area=6171148.7;  %Raster
% delta_axis=10;  %10 degrees

%%STARTING
cols=2160;
rows=900;
nondata=-32768;

Factor_Area=1/1000000;  %1000000000 m3 to BCM billion of cubic meters

folderMask='D:\02_DATA_D\WorldClim\';
load([folderMask,nameMask]);
eval(['MASK=',nameMask,';']);

% MASK=double(MASK);
% MASK(MASK==0)=NaN;

latmin=-60.0;
latmax=90.0;
longmin=-180.0;
longmax=180.0;

X_resolution=(longmax-longmin)/cols;
Y_resolution=(latmax-latmin)/rows;

% for period=pi:pf %loop
%     for path=ri:rf
% 
% period_opc=period; %period: 1)50, 2)70
% path_opc=path; %path: 1)26, 2)45, 3)60, 4)85

%loop

%%% (1.3)LOADING MATRICES

%%%STARTING
%Matrices-Models
model_mat={};
model_mat_differences={};
max_mat=[];    %1)States  2)Differences
min_mat=[];
sum_mat=[];
avg_mat=[];
count_mat=[];

max_mat2=[];    %1)States  2)Differences
min_mat2=[];
sum_mat2=[];
avg_mat2=[];
count_mat2=[];

I=find(MASK==1);
% Inan=find(isnan(MASK));
Inan=find(MASK==0);
land=sum(MASK(I));

%This because model (18) MPI-ESM-LR have values with errors
MASK=double(MASK);
MASKaux=ones(size(MASK));
MASKaux(691:end,:)=0;
MASKaux=MASK.*MASKaux;
Iaux=find(MASKaux==1);
landaux=sum(MASKaux(Iaux));
            
tic
 for v=variable_opc:variable_opc
    for p=period_opc:period_opc
       for r=path_opc:path_opc
                       
%Variable
if v==1
    val='pr';
    factor=1.0;
elseif v==2
    val='tm';
    factor=10.0;
elseif v==3
    val='tn';
    factor=10.0;
elseif v==4
    val='tx';
    factor=10.0;
elseif v==5
    val='bp';
    factor=1.0;
elseif v==6
    val='bt';
    factor=10.0;
end

%Period
if p==1
  period='50';
elseif p==2
  period='70';
end

%Path
 if r==1
   path='26';
elseif r==2
  path='45';
elseif r==3
  path='60';
elseif r==4
  path='85';
 end

 %Load matrix
 eval(['load ',char(39),folderbigMAT,suffix,val,period,path,'.mat',char(39),';']);
 
 for model_id=1:20
    for k=id_variable:id_variable
        opc=0;
    	if v<5
        	opc=1;
        else
            if v==5
                if k<=8
                	opc=1;
                end 
            elseif v==6
            	if k<=11
                	opc=1;    
                end
            end
        end
            
        if opc
            
            eval(['aux01(:,:)=',suffix,val,period,path,'(',int2str(model_id),',',int2str(k),',:,:);']);  
            
            
            %Temperatures from model 18) MPI-ESM-LR are wrong below
            %degree -30
            if and(r==1,and(model_id==18,and(v~=1,v~=5)))
                mask_mat=MASKaux;
                landval=landaux;
            else
                mask_mat=MASK;
                landval=land;
            end
            aux02=aux01;
            aux02(aux02==nondata)=0;
            aux02=double(aux02).*mask_mat/factor;
                      
            model_mat{ model_id,1}=aux02;
            
            max_mat(model_id,1)=max(aux02(:));
            if max_mat(model_id,1)==0
                min_mat(model_id,1)=0;
            else
                min_mat(model_id,1)=min(aux02(aux02~=0));
            end
            
            sum_mat(model_id,1)=sum(aux02(:));
%             count_mat(model_id,1)=sum(MASK(:));
            avg_mat(model_id,1)=sum_mat(model_id,1)/landval;
            
            
            %%%For differences
            if model_id==1
                aux03=mask_mat;
                model_mat_differences{model_id,1}=aux03;
            else
                if and(round(max_mat(model_id,1))==0,round(avg_mat(model_id,1))==0)
                    model_mat_differences{model_id,1}=zeros(size(aux03));    %No-data matrix
                else
                    aux03=100.*mask_mat.*(aux02-model_mat{1,1})./model_mat{1,1};
                    model_mat_differences{model_id,1}=aux03;
                end
            end
            
            max_mat2(model_id,1)=max(aux03(:));
            if max_mat(model_id,1)==0
               min_mat2(model_id,1)=0;
            else
               min_mat2(model_id,1)=min(aux03(aux03~=0));
            end
            
            sum_mat2(model_id,1)=sum(aux03(:));
%             count_mat2(model_id,1)=sum(MASK(:));
            avg_mat2(model_id,1)=sum_mat2(model_id,1)/land;
            
        end
    end
 end
 
 %Close matrix
 %eval(['save(',char(39),folderbigMAT,suffix,val,period,path,'.mat',char(39),',',char(39),folderbigMAT,suffix,val,period,path,char(39),');']);
 eval(['clear ',suffix,val,period,path,';']);
 
       end
    end
 end
toc
disp('End!!');

%%
%(1.4) PLOTTING 20 MAPS
%****Creating Axes*****
%**Limits
X_axis=zeros(cols+1,1); %size+zero
X_axis_default=zeros(cols+1,1);
Y_axis=zeros(rows+1,1);
Y_axis_default=zeros(rows+1,1);

%***Creates Label Axes (Totals): Middle***
X_axis(1,1)=longmin;
X_axis_default(1,1)=0;
for ii=2:(cols+1)
    X_axis(ii,1)=X_axis(ii-1,1)+X_resolution;
    X_axis_default(ii,1)=X_axis_default(ii-1,1)+1;
end
Y_axis(1,1)=latmin;
Y_axis_default(1,1)=0;
for jj=2:(rows+1)
    Y_axis(jj,1)=Y_axis(jj-1,1)+Y_resolution;
    Y_axis_default(jj,1)=Y_axis_default(jj-1,1)+1;
end
Y_axis=flip(Y_axis);  %flip axis labels

%X_labels
iii_count=0;
for iii=xmin:xmax+1
    if mod(round(X_axis(iii,1),1),delta_axis)==0       
      iii_count=iii_count+1;
    end
end
xxx_label=zeros(iii_count,1);
xxx=zeros(iii_count,1);

jjj_count=0;
for jjj=xmin:xmax+1
    if mod(round(X_axis(jjj,1),1),delta_axis)==0
      jjj_count=jjj_count+1;
      xxx_label(jjj_count,1)=round(X_axis(jjj,1),1);
      xxx(jjj_count,1)=jjj-xmin+1;
    end
end

%Y_labels
iii_count=0;
for iii=ymin:ymax+1
    if mod(round(Y_axis(iii,1),1),delta_axis)==0       
      iii_count=iii_count+1;
    end
end
yyy_label=zeros(iii_count,1);
yyy=zeros(iii_count,1);

jjj_count=0;
for jjj=ymin:ymax+1
    if mod(round(Y_axis(jjj,1),1),delta_axis)==0
      jjj_count=jjj_count+1;
      yyy_label(jjj_count,1)=round(Y_axis(jjj,1),1);
      yyy(jjj_count,1)=jjj-ymin+1;  
    end
end

%Subplot
iplot=4; %subplot rows
jplot=5; %subplot cols

%***STATES***

%Climits: limits of plotting data
if or(variable_opc==1,variable_opc==5)
    ClimXmin=0;  
else
    ClimXmin=min(min_mat(min_mat~=0));
end
[ClimXmax,Imax]=max(max_mat(:));
%ClimXmax=2000;
%ClimXmax=5000;

auxxx=9999999999999;
Imin=1;
for iiii=1:length(max_mat)
    if max_mat(iiii,1)~=0
    	if max_mat(iiii,1)<auxxx
        	auxxx=max_mat(iiii,1);
        	Imin=iiii;
        end
    end
end

%***Titles, subtitles & sizes***
%Title
str_supertitle=[variable_label{variable_opc,1},': ',subvariable_label{id_variable,variable_opc},', ',period_label{period_opc,1},', ',path_label{path_opc,1}];
size_supertitle=15;

%Subtitle
size_str_title=10;
%Xlabel
size_str_xlabel=8;

%Size labels
size_Axis_labels=6;

%Plotting till an value
%ClimXmax=800;

%20 models
%h(1)=figure(1);
figure;
set(gcf,'color','w');

for i=1:20
subplot(iplot,jplot,i);

%Ploting Matrices 
aux=model_mat{i,1};

%with imagesc
% imagesc(aux(ymin:ymax,xmin:xmax));

%with pcolor
aux(Inan)=nan;
pcolor(rot90(aux(ymin:ymax,xmin:xmax),-1)');
shading flat; %Without grid
axis square

%Some matrices do not have data
if and(round(max_mat(i,1))==0,round(avg_mat(i,1))==0)
    %with pcolor
    subplot(iplot,jplot,i);
    aux=zeros(size(rot90(aux(ymin:ymax,xmin:xmax),-1)'));
    aux(:)=nan;
    pcolor(aux);
    shading flat; %Without grid
    
    text(round(((xmax-xmin)/2)),round((ymax-ymin))/2,'No data');
    
    axis off;
else
    axis on;
    
    str_total=num2str(round(round(avg_mat(i,1)).*mask_Area*Factor_Area));
    str_max=num2str(round(max_mat(i,1)));
    str_avg=num2str(round(avg_mat(i,1)));
    str_min=num2str(round(min_mat(i,1)));
    
    if or(variable_opc==1,variable_opc==5)
        sublabel_units=' (mm)';    
        str_xlabel=['Total = ',str_total,' bcm, Avg = ',str_avg,sublabel_units];
    else
        sublabel_units=' (°C)';
        str_xlabel=['Max = ',str_max,', Avg = ',str_avg,', Min = ',str_min,sublabel_units];
    end
    
    str_xlabel_color='k';
    if i==Imax
        str_xlabel_color='r';
    elseif i==Imin
        str_xlabel_color='b';
    end
    xlabel(str_xlabel,'FontSize',size_str_xlabel,'Color',str_xlabel_color);
end

%limits of plotting data
set(gca, 'CLim', [ClimXmin, ClimXmax]);
set(gca,'XTick',xxx,'YTick',yyy);

%with imagesc
%set(gca,'XTickLabel',xxx_label,'YTickLabel',yyy_label,'FontSize',size_Axis_labels);

%with pcolor
set(gca,'XTickLabel',xxx_label,'YTickLabel',flipud(yyy_label),'FontSize',size_Axis_labels);

%Subtitles
str_title=[num2str(i),') ',model_label{i,1}];
title(str_title,'FontSize',size_str_title); 

end

%with imagesc
% %background and zeros: WHITE
% myColorMap = jet; % Make a copy of jet.
% % Assign white (all 1's) to black (the first row in myColorMap).
% myColorMap(1, :) = [1 1 1];
% colormap(myColorMap); % Apply the colormap

if or(variable_opc==1,variable_opc==5)
    colormap(flipud(jet(32)))
else
    colormap(jet(32));
end
%SuperTitle
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5,1,str_supertitle,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',size_supertitle);

hb = colorbar;
hb.Label.String = '';
set(hb,'Position',[0.92 .25 0.02 0.5],'Ticks',[0 1],'TickLabels',[round(ClimXmin) round(ClimXmax)]);

% end  %Loop
% end %loop


%%
%***PLOTTING DIFFERENCES
h(2)=figure(2);

%Climits: limits of plotting data
ClimXmin=min(min_mat2(min_mat2~=0));
[ClimXmax,Imax]=max(max_mat2(:));

auxxx=9999999999999;
Imin=1;
for iiii=1:length(max_mat2)
    if max_mat2(iiii,1)~=0
    	if max_mat2(iiii,1)<auxxx
        	auxxx=max_mat2(iiii,1);
        	Imin=iiii;
        end
    end
end

%***Titles, subtitles & sizes***
%Title
str_supertitle=[variable_label{variable_opc,1},': ',subvariable_label{id_variable,variable_opc},', ',period_label{period_opc,1},', ',path_label{path_opc,1}];
size_supertitle=15;

%Subtitle
size_str_title=10;
%Xlabel
size_str_xlabel=8;

%Size labels
size_Axis_labels=6;

%20 models
for i=1:20
subplot(iplot,jplot,i);

%Ploting Matrices 
aux=model_mat_differences{i,1};
imagesc(aux(ymin:ymax,xmin:xmax));

%Some matrices do not have data
if and(round(max_mat(i,1))==0,round(avg_mat(i,1))==0)
    
    text(round(((xmax-xmin)/2)),round((ymax-ymin))/2,'No data');
    axis off;
else
   
    axis on;
    str_total=num2str(round(round(avg_mat2(i,1)).*mask_Area*Factor_Area));
    str_max=num2str(round(max_mat2(i,1)));
    str_avg=num2str(round(avg_mat2(i,1)));
    str_min=num2str(round(min_mat2(i,1)));
    str_xlabel=['Total = ',str_total,', Max = ',str_max,', Avg = ',str_avg,', Min = ',str_min,''];
    %str_xlabel=['Max = ',str_max,', Min = ',str_min];
    str_xlabel_color='k';
    if i==Imax
        str_xlabel_color='r';
    elseif i==Imin
        str_xlabel_color='b';
    end
    xlabel(str_xlabel,'FontSize',size_str_xlabel,'Color',str_xlabel_color);
end

ClimXmin=-100;
ClimXmax=100;
%limits of plotting data
set(gca, 'CLim', [ClimXmin, ClimXmax]);
%
set(gca,'XTick',xxx,'YTick',yyy);
set(gca,'XTickLabel',xxx_label,'YTickLabel',yyy_label,'FontSize',size_Axis_labels);

%Subtitles
str_title=[num2str(i),') ',model_label{i,1}];
title(str_title,'FontSize',size_str_title); 

end

% %background and zeros: WHITE
% set(gcf,'color','w');
% myColorMap = jet; % Make a copy of jet.
% % Assign white (all 1's) to black (the first row in myColorMap).
% myColorMap(1, :) = [0 1 1];
% colormap(myColorMap); % Apply the colormap

% colormap(jet);

%SuperTitle
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5,1,str_supertitle,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',size_supertitle);

hb = colorbar;
hb.Label.String = '';
set(hb,'Position',[0.92 .25 0.02 0.5],'Ticks',[0 1],'TickLabels',[round(ClimXmin) round(ClimXmax)]);

%%
%***Save figure***
folderbigFIGURES='D:\02_DATA_D\WorldClim\CMIP5\FIGURES\';

%s:states
name_figure=[folderbigFIGURES,suffix,val,period,path,'_',num2str(id_variable),'_s.fig'];
savefig(h(1),name_figure)

%d:differences
name_figure=[folderbigFIGURES,suffix,val,period,path,'_',num2str(id_variable),'_d.fig'];
savefig(h(2),name_figure)

%%
%***********************
%(2.1) Choosing variable
variable_opc = 2; %val: 1)pr, 2)tm, 3)tn, 4)tx, 5)bp, 6)bt
pi=1; %period: 1)50, 2)70
pf=2;
ri=2; %path: 1)26, 2)45, 3)60, 4)85
rf=2;
vi=1;  %1-12) 
vf=12;

nondata=-32768;

%%%STARTING
%Matrices-Models
model_mat={};
max_mat=[]; 
min_mat=[];
sum_mat=[];
avg_mat=[];
sum_models_mat=[];
avg_models_mat=[];

I=find(MASK==1);
% Inan=find(isnan(MASK));
Inan=find(MASK==0);
land=sum(MASK(I));

%This because model (18) MPI-ESM-LR have values with errors
MASK=double(MASK);
MASKaux=ones(size(MASK));
MASKaux(691:end,:)=0;
MASKaux=MASK.*MASKaux;
Iaux=find(MASKaux==1);
landaux=sum(MASKaux(Iaux));

model_n=1;
tic
 for v=variable_opc:variable_opc
    for p=pi:pf
       for r=ri:rf
                       
%Variable
if v==1
    val='pr';
    factor=1.0;
elseif v==2
    val='tm';
    factor=10.0;
elseif v==3
    val='tn';
    factor=10.0;
elseif v==4
    val='tx';
    factor=10.0;
elseif v==5
    val='bp';
    factor=1.0;
elseif v==6
    val='bt';
    factor=10.0;
end

%Period
if p==1
  period='50';
elseif p==2
  period='70';
end

%Path
 if r==1
   path='26';
elseif r==2
  path='45';
elseif r==3
  path='60';
elseif r==4
  path='85';
 end

 %Load matrix
 eval(['load ',char(39),folderbigMAT,suffix,val,period,path,'.mat',char(39),';']);
 
 for model_id=1:20
    if model_id==1 
        for k=vi:vf
                opc=0;
                if v<5
                    opc=1;
                else
                    if v==5
                        if k<=8
                            opc=1;
                        end 
                    elseif v==6
                        if k<=11
                            opc=1;    
                        end
                    end
                end

                if opc
                   
                    eval(['aux01(:,:)=',suffix,val,period,path,'(',int2str(model_id),',',int2str(k),',:,:);']);  
                    
                    aux02=aux01;
                    aux02(aux02==nondata)=0;
                    aux02=double(aux02).*MASK/factor;

                   % model_mat{ model_id,1}=aux02;

                    sum_mat(k,1)=sum(aux02(:));

                    avg_mat(k,1)=sum(aux02(:))/land;

                end

        end          

    else
        for k=vi:vf
            opc=0;
            if v<5
                opc=1;
            else
                if v==5
                    if k<=8
                        opc=1;
                    end 
                elseif v==6
                    if k<=11
                        opc=1;    
                    end
                end
            end

            if opc
                
                eval(['aux01(:,:)=',suffix,val,period,path,'(',int2str(model_id),',',int2str(k),',:,:);']);  

                %Temperatures from model 18) MPI-ESM-LR are wrong below
                %degree -30
                if and(r==1,and(model_id==18,and(v~=1,v~=5)))
                    mask_mat=MASKaux;
                    landval=landaux;
                else
                    mask_mat=MASK;
                    landval=land;
                end
                aux02=aux01;
                aux02(aux02==nondata)=0;
                aux02=double(aux02).*mask_mat/factor;
                   
               % model_mat{ model_id,1}=aux02;
                
                sum_models_mat(k,model_n)=sum(aux02(:));
                
                if sum(aux02(:))==0
                    avg_models_mat(k,model_n)=nan;
                else
                    avg_models_mat(k,model_n)=sum(aux02(:))/landval;
                end

            end

        end
        model_n=model_n+1;
    end
 end
 
 %Close matrix
 %eval(['save(',char(39),folderbigMAT,suffix,val,period,path,'.mat',char(39),',',char(39),folderbigMAT,suffix,val,period,path,char(39),');']);
 eval(['clear ',suffix,val,period,path,';']);
 
       end
    end
 end
toc
disp('End!!');

%%
%(2.2) PLOTTING
if variable_opc==1
    sublabel_units=' (mm)';    
else
    sublabel_units=' (°C)';
end
size_supertitle=12;
xLabels = {'J', 'F', 'M', 'A','M','J','J','A','S','O','N','D'};

ymin=0;
ymax=30;

grey = [0.8,0.8,0.8];

% %%%Plotting 19 models
% %%%BOXPLOT
figure()
set(gcf,'color','w');

supertitle=[variable_label{variable_opc,1},sublabel_units,', 19 GCM'];
annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
boxplot(avg_models_mat(:,:)');

ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
%ylim([ymin ymax]);
xlim([1 12]);
set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
set(gca,'XTickLabel', xLabels,'Fontsize', 10);

%%%HIETOGRAM
figure()
set(gcf,'color','w');

annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
for ii=1:length(avg_models_mat)/2
    ensemble_50=plot(avg_models_mat(:,ii),'c');
    hold on
end
for ii=((length(avg_models_mat)/2)+1):length(avg_models_mat)
    ensemble_70=plot(avg_models_mat(:,ii),'y');
    hold on
end

mean_ensemble=plot(nanmean(avg_models_mat(:,1:end)'),'r','LineWidth',1.5);
current=plot(avg_mat(:,1),'b','LineWidth',1.5);
hold off;

ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
%ylim([ymin ymax]);
xlim([1 12]);
set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
set(gca,'XTickLabel', xLabels,'Fontsize', 10);

legend([current,mean_ensemble,ensemble_50,ensemble_70],'Mean historical (1950-2000)','Mean ensemble (2041-2080)','2050 projection','2070 projection','Location','southoutside','Orientation','horizontal');


%%%PLOTTING ACCUMULATED
if variable_opc==1
    %%%Plotting 19 models aveg
    %%%BOXPLOT
    figure()
    set(gcf,'color','w');

    supertitle=[variable_label{variable_opc,1},sublabel_units,', 19 GCM'];
    annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
    boxplot(cumsum(avg_models_mat(:,:))');

    ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
    %ylim([ymin ymax]);
    xlim([1 12]);
    set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
    set(gca,'XTickLabel', xLabels,'Fontsize', 10);

    %%%HIETOGRAM
    figure()
    set(gcf,'color','w');

    annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
    for ii=1:length(avg_models_mat)/2
        ensemble_50=plot(cumsum(avg_models_mat(:,ii)),'c');
        hold on
    end
    for ii=((length(avg_models_mat)/2)+1):length(avg_models_mat)
        ensemble_70=plot(cumsum(avg_models_mat(:,ii)),'y');
        hold on
    end

    mean_ensemble=plot(cumsum(nanmean(avg_models_mat(:,1:end)')),'r','LineWidth',1.5);
    current=plot(cumsum(avg_mat(:,1)),'b','LineWidth',1.5);
    hold off;

    ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
    %ylim([ymin ymax]);
    xlim([1 12]);
    set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
    set(gca,'XTickLabel', xLabels,'Fontsize', 10);

    legend([current,mean_ensemble,ensemble_50,ensemble_70],'Mean historical (1950-2000)','Mean ensemble (2041-2080)','2050 projection','2070 projection','Location','southoutside','Orientation','horizontal');
        
end


%%
%%%Plotting 19 models per RPC
%%BOXPLOT

k=0;
for j=1:2
figure()
set(gcf,'color','w');

supertitle=[variable_label{variable_opc,1},sublabel_units,', ',period_label{j,1},', 19 GCM'];
annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
     
    for i=1:4

        subplot(2,2,i)
        boxplot(avg_models_mat(:,k+(i-1)*19+1:k+i*19)');

        ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
        %ylim([ymin ymax]);
        
        set(gca,'XTickLabel', xLabels,'Fontsize', 10);
        title(path_label{i,1});
    end
    k=k+76;
end

%%%HIETOGRAM
k=0;
for j=1:2
figure()
set(gcf,'color','w');

supertitle=[variable_label{variable_opc,1},sublabel_units,', ',period_label{j,1},', 19 GCM'];
annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);

    for i=1:4
        subplot(2,2,i)
        for ii=k+(i-1)*19+1:k+i*19
            
            ensemble=plot(avg_models_mat(:,ii),'Color',grey);
            hold on

        end
        mean_ensemble=plot(nanmean(avg_models_mat(:,k+(i-1)*19+1:k+i*19)'),'r','LineWidth',1.5);
        current=plot(avg_mat(:,1),'b','LineWidth',1.5);
        hold off;
        
        ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
        %ylim([ymin ymax]);
        xlim([1 12]);
        set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
        set(gca,'XTickLabel', xLabels,'Fontsize', 10);
        title(path_label{i,1});
        
            
    end
    k=k+76;
end
%legend([current,mean_ensemble,ensemble],'Mean historical (1950-2000)','Mean ensemble','Projection','Location','southoutside','Orientation','horizontal');

%%%PLOTTING ACCUMULATED
if variable_opc==1    
    %%%Plotting 19 models per RPC
    %%BOXPLOT

    k=0;
    for j=1:2
    figure()
    set(gcf,'color','w');

    supertitle=[variable_label{variable_opc,1},sublabel_units,', ',period_label{j,1},', 19 GCM'];
    annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);
        
        for i=1:4

            subplot(2,2,i)
        
            boxplot(cumsum(avg_models_mat(:,k+(i-1)*19+1:k+i*19))');

            ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
            %ylim([ymin ymax]);

            set(gca,'XTickLabel', xLabels,'Fontsize', 10);
            title(path_label{i,1});
        end
        k=k+76;
    end


    %%%HIETOGRAM
    k=0;
    for j=1:2
    figure()
    set(gcf,'color','w');

    supertitle=[variable_label{variable_opc,1},sublabel_units,', ',period_label{j,1},', 19 GCM'];
    annotation('textbox', [0 0.9 1 0.1],'String', supertitle,'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',size_supertitle);

        for i=1:4

            subplot(2,2,i)
            
            for ii=k+(i-1)*19+1:k+i*19

                ensemble=plot(cumsum(avg_models_mat(:,ii)),'Color',grey);
                hold on

            end
            mean_ensemble=plot(cumsum(nanmean(avg_models_mat(:,k+(i-1)*19+1:k+i*19)')),'r','LineWidth',1.5);
            current=plot(cumsum(avg_mat(:,1)),'b','LineWidth',1.5);
            hold off;

            ylabel([variable_label{variable_opc,1},sublabel_units], 'Fontsize', 10)
            %ylim([ymin ymax]);
            xlim([1 12]);
            set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12]);
            set(gca,'XTickLabel', xLabels,'Fontsize', 10);
            title(path_label{i,1});

        end
        k=k+76;
    end
    %legend([current,mean_ensemble,ensemble],'Mean historical (1950-2000)','Mean ensemble','Projection','Location','southoutside','Orientation','horizontal');
    
    
end


