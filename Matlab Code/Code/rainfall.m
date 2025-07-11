clear all
close all
clc

threshold = 2.5;                                                             % rainfall threshold







%% Loading netCDF file into Matlab Workspace
path = 'C:\Users\gess-icwar\OneDrive - Indian Institute of Science\Desktop\Project\Matlab Code\Rainfall Files\_Clim_Pred_LRF_New_RF25_IMD0p252019.nc';
file = ncinfo(path);
file.Variables.Name;
ncID = netcdf.open(path);

rainid = netcdf.inqVarID(ncID,"RAINFALL");
lonid = netcdf.inqVarID(ncID,"LONGITUDE");
latid = netcdf.inqVarID(ncID,"LATITUDE");

% Variavbles

rain = netcdf.getVar(ncID,rainid);
lat = netcdf.getVar(ncID,latid);
lon = netcdf.getVar(ncID,lonid);
 


index =1;
for p = 1:135
    for q = 1:129
%% Getting variables of interest


i = 1; wd = 0; drdy =0; spell = 0; drspell = 0; j = 1;                     %% Grid = p,q 
a = squeeze(rain(p,q,:));                                                  %% p = Longtitude ; q = latitude
if a(:,1) ~= -999
wd_data1 = zeros(365,3);
drdy_s_len = zeros(365,1);
for k = 1:365   
    if a(k,1) > threshold
        wd = wd+1;
        if k == 1
            spell = 1;
            wd_data1(i,1) = wd_data1(i,1) + a(k,1);
            wd_data1(i,2) = spell;
            wd_data1(i,3) = 1;
             i =  i + 1;
        else
            if a(k-1,1) > threshold
                spell = spell +1;
                wd_data1(i-1,1) = wd_data1(i-1,1) + a(k,1);
                wd_data1(i-1,2) = spell;
                wd_data1(i-1,3) = k;
                 
            else
                spell = 1;
                wd_data1(i,1) = wd_data1(i,1) + a(k,1);
                wd_data1(i,2) = 1;
                wd_data1(i,3) = k;
                i = i+1;
            end 
        end
    else 
        drdy = drdy + 1;
        if k == 1
            drspell = drspell + 1;
            drdy_s_len(j,1) = drspell;
            drdy_s_len(j,2) = 1;
            j = j + 1;
        else
            if a(k-1,1) < threshold
                if j == 1
                    drspell = drspell + 1;
                    drdy_s_len(j,1) = drspell;
                    drdy_s_len(j,2) = k;
                else
                     drspell = drspell + 1;
                     drdy_s_len(j-1,1) = drspell;
                     drdy_s_len(j-1,2) = k;
                end
            else
                if j == 1 && nnz(wd_data1(:,2)) < threshold
                    drspell = 1; 
                    drdy_s_len(j,1) = drspell;
                    drdy_s_len(j,2) = k;
                    
                else
                    drspell = 1; 
                    drdy_s_len(j+1,1) = drspell;
                    drdy_s_len(j+1,2) = k;
                    j = j + 1;
                end
            end
         end
    end
end

if drdy == 365
    wd_data1(:,1) = 0;
    wd_data1(:,2) = 0;
    wd_data1(:,3) = 0;
end

%% Corrections to days data
wd_data1(:,4) = wd_data1(:,3) - wd_data1(:,2) + 1;
drdy_s_len(:,3) = drdy_s_len(:,2) - drdy_s_len(:,1) + 1;

%% Removing end zeros
if wd_data1(1,1)  == 0 
    wd_data(1,:,index) = [0 0 0 0 lon(p) lat(q)];
else
    x = 1;
    while wd_data1(x,1) ~= 0
            wd_data(x,:,index) = [wd_data1(x,:) lon(p) lat(q)];
            x = x+1;
    end
end

if drdy_s_len(1,1) == 365
    dry_data(1,:,index) = [365 365 1 lon(p) lat(q)];
else
    y = 1;
    while drdy_s_len(y,1) ~= 0
            dry_data(y,:,index) = [drdy_s_len(y,:) lon(p) lat(q)];
            y = y+1;
    end
end

    %% Verifying sum of rainfall & no of wet & dry days
r_cum = 0;
for i =1:365
    if a(i,1) > threshold
        r_cum = r_cum + a(i,1);
    end 
end
verify(index,1) = sum(wd_data(:,1,index)) - r_cum ;
verify(index,2) = sum(wd_data(:,2,index)) - wd ;
verify_dry(index,1) = sum(drdy_s_len(:,1)) - drdy;
index = index + 1;
end    
end
    
end
