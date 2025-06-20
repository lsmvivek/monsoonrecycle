clear evp
sinificance_level = 0.2;
start_date = 150;
end_date = 288;








%%   Loading Evaporation netCDF file ( clipped-regridded using Xarray Package in Python ) to Matlab Workspace 
path = 'C:\Users\gess-icwar\OneDrive - Indian Institute of Science\Desktop\Project\Matlab Code\Evaporation Files\Clipped Evp\saved_on_disk_2019.nc';
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

%% Getting exact cumulative evaporation correponding to dry spells

  %  for k =1:length(dry_data(:,1,17415))


           
for k = 1:4964
    for i = 1:17415
        if (dry_data(1,4,k) == evapo(1,2,i) && dry_data(1,5,k) == evapo(1,3,i)  )
            index3 = 1;
            for j = 1:length(dry_data(:,1,k))    
                if (dry_data(j,1,k) ~= 0 && dry_data(j,3,k) > start_date && dry_data(j,3,k)< end_date )     % May 31 till October 30
                    evp(index3,2,k)= wd_data(j,1,k);                              % X1 wet spell cumulative
                    evp(index3,3,k)= wd_data(j,2,k);                              % X2 wet spell length
                    evp(index3,1,k) = sum(evapo(dry_data(j,3,k):dry_data(j,2,k),1,i));   % Evaporation in following dry spell
                    evp(index3,4,k) = dry_data(1,4,k);                           % longitude
                    evp(index3,5,k) = dry_data(1,5,k);                           % latitude
                    evp(index3,6,k) = wd_data(j,4,k);                            % start of WS
                    evp(index3,7,k) = wd_data(j,3,k);                            % end of WS
                    evp(index3,8,k) = dry_data(j,1,k);                           % length of DS
                    evp(index3,9,k) = dry_data(j,2,k);                           % end of DS
                    evp(index3,10,k) = dry_data(j,3,k);                          % start of DS
                    index3 = index3 + 1;
                end
            end
        end    
    end
end

%% Averaging the data
for i = 1:length(evp(1,1,:))
    alpha(:,1,i) = evp(:,1,i);
    alpha(:,2,i) = evp(:,2,i);
    alpha(:,3,i) = evp(:,3,i);
end



%% Calculation of correlation 

corrmat = zeros(length(alpha(1,1,:)),2);
pvalmat = zeros(length(alpha(1,1,:)),2);
for i = 1:length(evp(1,1,:))
    if  evp(1,1,i) ~= 0 && dry_data(1,2,i) ~= 365        % NaN values of evaporation & full dry year grid excluded              
        y1 = nonzeros(alpha(1:nnz(evp(:,1,i)),1,i));
        x1 = nonzeros(alpha(1:nnz(evp(:,2,i)),2,i));
        x2 = nonzeros(alpha(1:nnz(evp(:,3,i)),3,i));
        if length(y1)==length(x1)
            if evp(1,10,i)== start_date + 1
                if length(y1) ==1
                    [corrmat(i,1) pvalmat(i,1)] = corr(y1,x1,"tail","both","type","Kendall");
                    [corrmat(i,2) pvalmat(i,2)] = corr(y1,x2,"tail","both","type","Kendall");
                else
                    [corrmat(i,1) pvalmat(i,1)] = corr(y1(2:end),x1(1:end-1),"tail","both","type","Kendall");
                    [corrmat(i,2) pvalmat(i,2)] = corr(y1(2:end),x2(1:end-1),"tail","both","type","Kendall");    
                end
            else
                    [corrmat(i,1) pvalmat(i,1)] = corr(y1,x1,"tail","both","type","Kendall");
                    [corrmat(i,2) pvalmat(i,2)] = corr(y1,x2,"tail","both","type","Kendall");
            end
        elseif length(y1) == length(x1) -1   
               [corrmat(i,1) pvalmat(i,1)] = corr(y1,x1(1:end-1),"tail","both","type","Kendall");
               [corrmat(i,2) pvalmat(i,2)] = corr(y1,x2(1:end-1),"tail","both","type","Kendall");
        elseif length(y1) -1 == length(x1) 
            if length(y1) == 1
               [corrmat(i,1) pvalmat(i,1)] = corr(y1(1,1),x1,"tail","both","type","Kendall");
               [corrmat(i,2) pvalmat(i,2)] = corr(y1(1,2),x2,"tail","both","type","Kendall");
            else
               [corrmat(i,1) pvalmat(i,1)] = corr(y1(2:end),x1,"tail","both","type","Kendall");
               [corrmat(i,2) pvalmat(i,2)] = corr(y1(2:end),x2,"tail","both","type","Kendall");
            end
        end
    else
         corrmat(i,1) = -999;
         pvalmat(i,1) = -999;
         corrmat(i,2) = -999;
         pvalmat(i,2) = -999;
    end
end

for i = 1:length(alpha(1,1,:))
    N(i) = nnz(alpha(:,1,i));
    signi(i) = 1.96/sqrt(N(i) - 3);
end


%% Removing Nan results, at grids if anyone correlation was not available it was ignored


index4 = 1;
misseding = 0;
for i = 1:length(evp(1,1,:))
    if  corrmat(i,1) ~= -999 && ~isnan(corrmat(i,1)) && corrmat(i,2) ~= -999 && ~isnan(corrmat(i,2))
        all_vals(index4,1:2) = corrmat(i,:);                                
        all_p(index4,1:2)=pvalmat(i,:);
        all_vals(index4,3) = evp(1,4,i);                                    % Longitude
        all_vals(index4,4) = evp(1,5,i);                                    % Latitude 
        index4 = index4 + 1;
        if  abs(corrmat(i,1)) == 1 && abs(corrmat(i,2)) == 1
            all_p(index4-1,1:2)= [0 0];
        end
    else
        misseding = misseding + 1;
    end
end
% %% Figures     % first two are of Rajarshi sir's
% count_sign_x1 = 0;
% count_sign_x2 = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RDJ X1
% figure;
% plot(all_vals(:,1))
% yline(0.4)
% yline(-0.4)
% hold on
% for i = 1:length(all_vals)
%     if abs(all_vals(i,1)) > signi(i) 
%         count_sign_x1 = count_sign_x1 + 1;
%          plot(i,all_vals(i,1),"*r");
%     else
%         plot(i,all_vals(i,1));
%     end
% 
% end
% 
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RDB X2
% figure;
% plot(all_vals(:,2))
% yline(0.4)
% yline(-0.4)
% hold on
% for i = 1:length(all_vals)
%     if abs(all_vals(i,2)) > signi(i) 
%         count_sign_x2 = count_sign_x2 + 1;
%          plot(i,all_vals(i,2),"*r");
%     else
%         plot(i,all_vals(i,2));
%     end
% end
% 
% hold off
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pval plots X1
figure;
plot(all_vals(:,1));
yline(0.4)
yline(-0.4)

hold on 
count_X1=0;
for i=1:length(all_p)
    if all_p(i,1)<=sinificance_level
       plot(i,all_vals(i,1),"*r");
       count_X1=count_X1+1;
    else
        plot(i,all_vals(i,1));
    end
end
hold off




%%%%%%%%%%%%%%%%%%%%%%%%
% pval plot X2
figure;
plot(all_vals(:,2));
yline(0.4)
yline(-0.4)

hold on 
count_X2=0;

for i=1:length(all_p)
        if all_p(i,2)<=sinificance_level
           plot(i,all_vals(i,2),"*r");
           count_X2=count_X2+1;
        else
            plot(i,all_vals(i,2));
        end
end
hold off

%% Significant only

% Change here for X1 & X2 both or separate

index5 = 1;
for j=1:length(all_p(:,1))
    for i = 1:length(evp(1,1,:))
       if all_p(j,1)<= sinificance_level  && all_vals(j,3) == evp(1,4,i) && all_vals(j,4) == evp(1,5,i)

           beta_one(index5,1) = evp(1,4,i);                                   % Longitude
           beta_one(index5,2) = evp(1,5,i);                                   % Latitude
           beta_one(index5,3) = all_vals(j,1);                                % X1 corr coefficeints       
           index5 = index5 + 1;
       end
    end
end



index6 = 1;
for j=1:length(all_p(:,2))
    for i = 1:length(evp(1,1,:))
       if all_p(j,2)<=sinificance_level  && all_vals(j,3) == evp(1,4,i) && all_vals(j,4) == evp(1,5,i)

           beta_two(index6,1) = evp(1,4,i);                                   % Longitude
           beta_two(index6,2) = evp(1,5,i);                                   % Latitude    
           beta_two(index6,3) = all_vals(j,2);                                % X2 corr coefficients
           index6 = index6 + 1;
       end
    end
end
%% Take beta for significant results with lat & lon
 % Take all_vals for all results with lat & lon




