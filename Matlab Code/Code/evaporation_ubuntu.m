clear evp
%%   Loading Evaporation netCDF file ( clipped-regridded using Xarray Package in Python ) to Matlab Workspace 
path = '/home/gess-icwar/Desktop/Evaporation Files/python code/saved_on_disk_2011.nc';
file = ncinfo(path);
file.Variables.Name;

ncID = netcdf.open(path);

eid = netcdf.inqVarID(ncID,"E");
dlonid = netcdf.inqVarID(ncID,"lon");
dlatid = netcdf.inqVarID(ncID,"lat");
bandid = netcdf.inqVarID(ncID,"time");
 
e = netcdf.getVar(ncID,eid);
dlat = netcdf.getVar(ncID,dlatid);
dlon = netcdf.getVar(ncID,dlonid);
time = netcdf.getVar(ncID,bandid);


%% Getting exact cumulative evaporation correponding to dry spells
% Saved daily evaporation gridwise
index2 = 1;
for p = 1:135
    for q = 1:129                                                                                                                                
        b = squeeze(e(p,q,:));
        evapo(:,1,index2) = b;
        evapo(:,2,index2) = dlon(p);
        evapo(:,3,index2) = dlat(q);
        index2 = index2 + 1;
    end
end

% Getting exact cumulative evaporation correponding to dry spells

  %  for k =1:length(dry_data(:,1,17415))


           
for k = 1:4964
    for i = 1:17415
        if (dry_data(1,4,k) == evapo(1,2,i) && dry_data(1,5,k) == evapo(1,3,i))
            index3 = 1;
            for j = 1: length(wd_data(:,1,1))    
                if (dry_data(j,1,k) ~= 0 && dry_data(j,3,k) > 180 && dry_data(j,3,k)<300)
                    evp(index3,2,k)= wd_data(j,1,k);                              % X1 wet spell length
                    evp(index3,3,k)= wd_data(j,2,k);                              % X2 wet spell cumulative
                    evp(index3,1,k) = sum(evapo(dry_data(j,3,k):dry_data(j,2,k),1,i));   % Evaporation in following dry spell
                    evp(index3,4,k) = dry_data(1,4,k); 
                    evp(index3,5,k) = dry_data(1,5,k);
                    evp(index3,6,k) = wd_data(j,4,k);
                    evp(index3,7,k) = wd_data(j,3,k);
                    evp(index3,8,k) = dry_data(j,1,k);
                    evp(index3,9,k) = dry_data(j,2,k);
                    evp(index3,10,k) = dry_data(j,3,k);
                    index3 = index3 + 1;
                end
            end
        end    
    end
end
     
corrmat = zeros(4964,2);
pvalmat = zeros(4964,2);
for i = 1:4964
    if ~isnan(evp(1,1,i)) ~= 0 && dry_data(1,2,i) ~= 365
        y1 = nonzeros(evp(:,1,i));
        x1 = nonzeros(evp(:,2,i));
        x2 = nonzeros(evp(:,3,i));
        if length(y1)==length(x1)
            [corrmat(i,1) pvalmat(i,1)] = corr(y1,x1,"tail","both","type","Spearman");
            [corrmat(i,2) pvalmat(i,2)] = corr(y1,x2,"tail","both","type","Spearman");
        elseif length(y1) == length(x1) -1   
               [corrmat(i,1) pvalmat(i,1)] = corr(y1,x1(2:end),"tail","both","type","Spearman");
               [corrmat(i,2) pvalmat(i,2)]  = corr(y1,x2(2:end),"tail","both","type","Spearman");
        elseif length(y1) -1 == length(x1)   
               [corrmat(i,1) pvalmat(i,1)] = corr(y1(2:end),x1,"tail","both","type","Spearman");
               [corrmat(i,2) pvalmat(i,2)] = corr(y1(2:end),x2,"tail","both","type","Spearman");
        end
    else
         corrmat(i,1) = -999;
         pvalmat(i,1) = -999;
         corrmat(i,2) = -999;
         pvalmat(i,2) = -999;
    end
end

negatives = 0;
positives = 0;
missing = 0;

negsum_a = 0;
negsum_b = 0;
possum_a = 0;
possum_b = 0;

index4 = 1;

for i = 1:4964
    if corrmat(i,1) < 0 && corrmat(i,1) > -1
        negatives = negatives + 1;
        negsum_a = negsum_a + corrmat(i,1);
        negsum_b = negsum_b + corrmat(i,2);
        all_vals(index4,1:2) = corrmat(i,:);
        all_p(index4,:)=pvalmat(i,:);
        all_vals(index4,3) = evp(1,4,i);
        all_vals(index4,4) = evp(1,5,i);
        index4 = index4 + 1;
    elseif corrmat(i,1) > 0  
        all_vals(index4,1:2) = corrmat(i,:);
        all_p(index4,:)=pvalmat(i,:);
        all_vals(index4,3) = evp(1,4,i);
        all_vals(index4,4) = evp(1,5,i);
        positives = positives + 1;
        possum_a = possum_a + corrmat(i,1);
        possum_b = possum_b + corrmat(i,2);
        index4 = index4 + 1;
    else
        missing = missing + 1;
    end
end
figure;
plot(all_vals(:,1))
yline(0.2,'r')
yline(-0.2,'r')
hold on 
count_p = 0;
for i=1:length(all_vals)
if all_p(i,1)<=0.2
    plot(i,all_vals(i,1),"*k")
    count_p = count_p + 1;
else
    plot(i,0)
end
end
hold off


figure;
plot(all_vals(:,2))
yline(0.2,'r')
yline(-0.2,'r')
for l=1:4751
data(l,:) = [all_vals(l,1:4), all_p(l,1:2)];
end
