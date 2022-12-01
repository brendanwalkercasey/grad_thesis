   
%=========================================================================
% Description: The code below performs correlation analysis on Palmer and Standard Precipitiation Drought
% Indicies and Wildfire Total Acres Burned (1984-2019).
  
% Data:
%       Wildfire - https://www.mtbs.gov/product-descriptions
%       Drought - https://www.ncei.noaa.gov/pub/data/cirs/climdiv/drought-readme.txt

% Input:
%       dataRead.m - read drought and wildfire variables
%       PDSI - Palmer Drought Severity Index
%       PHDI - Palmer Hydrologic Drought Index
%       PMDI - Palmer Modified Drought Index
%       ZNDX - Palmer-Z Index
%       SP01 - Standard Precipitation Index
%       ACRES - Total Burned Acres (perim. delin.)

% Output: 
%       Kendall's tau, Spearman Rank, and Pearson's correlation coefficient
%  
% Author: Brendan Casey
% Assign: Graduate Thesis (M.S. in Geographic and Cartographic Sciences)
% Defense Date: 20220818
% Modification Records: 20221130 (cleaned code, added descrip.)
%                       20221201 (cleaned code, rev. descrip.)
%=========================================================================

%% Steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert all arrays to Non-Zero Arrays (for better analysis on drought w/ wildfire events)
% source - (https://www.mathworks.com/matlabcentral/answers/60877-remove-rows-from-a-matrix-on-a-specific-condition)
    % Enter these sequentially:
    % 1)
    AllData = [acres, pdsi_all, phdi_all, pmdi_all, zndx_all, sp01_all];
    % 2)
    noZero_Data = AllData;
    % 3)
    noZero_Data(noZero_Data(:, 1)== 0, :)= [];
    % 4)
    writematrix(noZero_Data,'all_data_noZeroWildfire.txt');
    % 5)
    fireNoZ = load('all_data_noZeroWildfire.txt');
    % 6)
%     fire_NoZ = fireNoZ(:,1);
%     pdsi_all_noZ = fireNoZ(:,2);
%     phdi_all_noZ = fireNoZ(:,3);
%     pmdi_all_noZ = fireNoZ(:,4);
%     zndx_all_noZ = fireNoZ(:,5);
%     sp01_all_noZ = fireNoZ(:,6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Prob plots:
% t = tiledlayout(2,3);
% %figure()
% nexttile
% normplot(phdi_all)
% title("PHDI")
% %figure()
% nexttile
% normplot(pdsi_all)
% title("PDSI")
% %figure()
% nexttile
% normplot(pmdi_all)
% title("PMDI")
% %figure()
% nexttile
% normplot(zndx_all)
% title("ZNDX")
% %figure()
% nexttile
% normplot(sp01_all)
% title("SP01")
% %figure()
% nexttile
% normplot(acres)
% title("Burn Severity (Acres)")


%% Separate PHDI distributions (both normal distirbutions, but when together, they are binomial) (separated by positive and negative)
% phdi_pos = load('PHDI_PosVal.txt');
% phdi_neg = load('PHDI_NegVal.txt');

%     noZero_Data = phdi_pos;
%     noZero_Data(noZero_Data(:, 3)== 0, :)= [];
%     writematrix(noZero_Data,'phdi_pos_noZeroWildfire.txt');
%     phdiPosfireNoZ = load('phdi_pos_noZeroWildfire.txt');

%     noZero_Data = phdi_neg;
%     noZero_Data(noZero_Data(:, 3)== 0, :)= [];
%     writematrix(noZero_Data,'phdi_neg_noZeroWildfire.txt');
%     phdiNegfireNoZ = load('phdi_neg_noZeroWildfire.txt');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Stats on diff. sets of variables:

    %% All Data (options)
%         allData = [fire_NoZ,pdsi_all_noZ];
        allData = [acres, pdsi_all];
%         allData = [acres, phdi_all];
%         allData = [acres, pmdi_all];
%         allData = [acres, zndx_all];
%         allData = [acres, sp01_all];


    %% Target Years
%         2014 - 2018 (3 consec. months top 25th percentile)
%             allData = [fireNoZ(218:256,1), fireNoZ(218:256,2)]; %psdi
%             allData = [fireNoZ(218:256,1), fireNoZ(218:256,3)]; %phdi
%             allData = [fireNoZ(218:256,1), fireNoZ(218:256,4)]; %pmdi
%             allData = [fireNoZ(218:256,1), fireNoZ(218:256,5)]; %zndx
%             allData = [fireNoZ(218:256,1), fireNoZ(218:256,6)]; %sp01
  
%         1987 (atleast one month in top 25th percentile above 600000)
%             allData = [fireNoZ(25:33,1), fireNoZ(25:33,2)]; %psdi
%             allData = [fireNoZ(25:33,1), fireNoZ(25:33,3)]; %phdi
%             allData = [fireNoZ(25:33,1), fireNoZ(25:33,4)]; %pmdi
%             allData = [fireNoZ(25:33,1), fireNoZ(25:33,5)]; %zndx
%             allData = [fireNoZ(25:33,1), fireNoZ(25:33,6)]; %sp01

%         2002-2003 (atleast one month in top 25th percentile above 600000)
%             allData = [fireNoZ(134:147,1), fireNoZ(134:147,2)]; %psdi
%             allData = [fireNoZ(134:147,1), fireNoZ(134:147,3)]; %phdi
%             allData = [fireNoZ(134:147,1), fireNoZ(134:147,4)]; %pmdi
%             allData = [fireNoZ(134:147,1), fireNoZ(134:147,5)]; %zndx
%             allData = [fireNoZ(134:147,1), fireNoZ(134:147,6)]; %sp01

%         2008 (atleast one month in top 25th percentile above 600000)
%             allData = [fireNoZ(181:187,1), fireNoZ(181:187,2)]; %psdi
%             allData = [fireNoZ(181:187,1), fireNoZ(181:187,3)]; %phdi
%             allData = [fireNoZ(181:187,1), fireNoZ(181:187,4)]; %pmdi
%             allData = [fireNoZ(181:187,1), fireNoZ(181:187,5)]; %zndx
%             allData = [fireNoZ(181:187,1), fireNoZ(181:187,6)]; %sp01

    %% Seasons (all years)
        %% Spring
%             allData = [tot_spr(:), pdsi_spr];
%             allData = [tot_spr(:), phdi_spr];
%             allData = [tot_spr(:), pmdi_spr];
%             allData = [tot_spr(:), zndx_spr];
%             allData = [tot_spr(:), sp01_spr];
            
            % 3 month avg
%                     Avg_tot_spr = mean(reshape(tot_spr, 3, []));
%                     Avg_pdsi_spr = mean(reshape(pdsi_spr, 3, []));
%                     Avg_pmdi_spr = mean(reshape(pmdi_spr, 3, []));
%                     Avg_phdi_spr = mean(reshape(phdi_spr, 3, []));
%                     Avg_zndx_spr = mean(reshape(zndx_spr, 3, []));
%                     Avg_sp01_spr = mean(reshape(sp01_spr, 3, []));
%                 allData = [Avg_tot_spr(:), Avg_sp01_spr(:)];
                
        %% Summer
%             allData = [tot_summ(:), pdsi_summ];
%             allData = [tot_summ(:), phdi_summ];
%             allData = [tot_summ(:), pmdi_summ];
%             allData = [tot_summ(:), zndx_summ];
%             allData = [tot_summ(:), sp01_summ];

            % 3 month avg
%                     Avg_tot_summ = mean(reshape(tot_summ, 3, []));
%                     Avg_pdsi_summ = mean(reshape(pdsi_summ, 3, []));
%                     Avg_pmdi_summ = mean(reshape(pmdi_summ, 3, []));
%                     Avg_phdi_summ = mean(reshape(phdi_summ, 3, []));
%                     Avg_zndx_summ = mean(reshape(zndx_summ, 3, []));
%                     Avg_sp01_summ = mean(reshape(sp01_summ, 3, []));
%                 allData = [Avg_tot_summ(:), Avg_sp01_summ(:)];
                
        %% Fall
%             allData = [tot_fall(:), pdsi_fall];
%             allData = [tot_fall(:), phdi_fall];
%             allData = [tot_fall(:), pmdi_fall];
%             allData = [tot_fall(:), zndx_fall];
%             allData = [tot_fall(:), sp01_fall];

            % 3 month avg
%                     Avg_tot_fall = mean(reshape(tot_fall, 3, []));
%                     Avg_pdsi_fall = mean(reshape(pdsi_fall, 3, []));
%                     Avg_pmdi_fall = mean(reshape(pmdi_fall, 3, []));
%                     Avg_phdi_fall = mean(reshape(phdi_fall, 3, []));
%                     Avg_zndx_fall = mean(reshape(zndx_fall, 3, []));
%                     Avg_sp01_fall = mean(reshape(sp01_fall, 3, []));
%                 allData = [Avg_tot_fall(:), Avg_sp01_fall(:)];
                
        %% Winter
%             allData = [tot_wint(:), pdsi_wint];
%             allData = [tot_wint(:), phdi_wint];
%             allData = [tot_wint(:), pmdi_wint];
%             allData = [tot_wint(:), zndx_wint];
%             allData = [tot_wint(:), sp01_wint];

            % 3 month avg
%                     Avg_tot_wint = mean(reshape(tot_wint, 3, []));
%                     Avg_pdsi_wint = mean(reshape(pdsi_wint, 3, []));
%                     Avg_pmdi_wint = mean(reshape(pmdi_wint, 3, []));
%                     Avg_phdi_wint = mean(reshape(phdi_wint, 3, []));
%                     Avg_zndx_wint = mean(reshape(zndx_wint, 3, []));
%                     Avg_sp01_wint = mean(reshape(sp01_wint, 3, []));
%                 allData = [Avg_tot_wint(:), Avg_sp01_wint(:)];
% %         
         %% Offset
%             allData = [Avg_tot_spr(:), Avg_sp01_summ(:)];
%             allData = [Avg_tot_summ(:), Avg_sp01_fall(:)];
%             allData = [Avg_tot_fall(:), Avg_sp01_wint(:)];


         %% Summer vs Fall (all years)    
            %allData = [tot_summ(:), pdsi_fall];
    %% Selected Seasons (high burn years)    
        %allData = [sixMon85_Fr, sixmon85_Dr];
         %% 1985 only
        %allData = [fire85, pdsi85(:)];
         %% 1987 only
        %allData = [fire87, pdsi87(:)];
    %% 1984 only
        %allData = [fire84, pdsi84(:)];
    %% 1984 (high burn severity only)
        %allData = [hi84, pdsi84(:)];
    %% 1985 (high burn severity only)
        %allData = [hi85, pdsi85(:)];
    %% El Nino Years    
        %allData = [elNinoFrYrs, elNinoPDSI(:)]; %pdsi
    %% El Nino Offset years (ElNino Year Drought vs Proceeding Year Fire)    
        %allData = [elNinoFrOff, elNinoPDSI(:)];
    %% Num. of Fire Events vs drought    
        %allData = [event, pdsi_all]; %pdsi
        %allData = [event, phdi_all]; %phdi
        %allData = [event, pmdi_all]; %pmdi
        %allData = [event, zndx_all]; %zndx
        %allData = [event, sp01_all]; %sp01
   %% Offset
        % all acres offset (prev year drought vs current year fire)
        %allData = [acres, pdsi_off1Yr];
        % one month (prev mon drought, curr mon fire)
%             allData = [acres, pdsi_1monOffset];
%             allData = [acres, phdi_1monOffset];
%             allData = [acres, pmdi_1monOffset];
%             allData = [acres, zndx_1monOffset];
%             allData = [acres, sp01_1monOffset];
        % two month
%             allData = [acres, pdsi_2monOffset];
%             allData = [acres, phdi_2monOffset];
%             allData = [acres, pmdi_2monOffset];
%             allData = [acres, zndx_2monOffset];
%             allData = [acres, sp01_2monOffset];
        % three month
%             allData = [acres, pdsi_3monOffset];
%             allData = [acres, phdi_3monOffset];
%             allData = [acres, pmdi_3monOffset];
%             allData = [acres, zndx_3monOffset];
%             allData = [acres, sp01_3monOffset];
    
    %% Lag (one year)
%             allData = [acres,pdsi_1yrOffset];
%             allData = [acres,phdi_1yrOffset];
%             allData = [acres,pmdi_1yrOffset];
%             allData = [acres,zndx_1yrOffset];
%             allData = [acres,sp01_1yrOffset];
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stats
%% Correlation tests (Rho only)
%     r_pearson = corr(allData,'type','Pearson')
%     r_spearman = corr(allData,'Type','Spearman')
%     r_kendall = corr(allData,'Type','Kendall')
% 

%% Correlation tests (P-Val and Rho)
[Prho, Ppval] = corr(allData, 'type', 'Pearson')
[Srho, Spval] = corr(allData, 'type', 'Spearman')
[Krho, Kpval] = corr(allData, 'type', 'Kendall')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra   
    
%{     
    
%% Multiple Regression
  
    %% Drought Indices
       
        % Total Acres Burned
            Y = acres;
            x1 = pdsi_all;
            x2 = phdi_all;
            x3 = pmdi_all;
            x4 = zndx_all;
            x5 = sp01_all;
            X = [x1, x2, x3, x4, x5];
             
        % PDSI  
            Y = pdsi_all;
            x1 = acres;
            x2 = phdi_all;
            x3 = pmdi_all;
            x4 = zndx_all;
            x5 = sp01_all;
            X = [x1, x2, x3, x4, x5];
            
        % PHDI  
            Y = phdi_all;
            x1 = pdsi_all;
            x2 = acres;
            x3 = pmdi_all;
            x4 = zndx_all;
            x5 = sp01_all;
            X = [x1, x2, x3, x4, x5];

        % PMDI  
            Y = pmdi_all;
            x1 = pdsi_all;
            x2 = phdi_all;
            x3 = acres;
            x4 = zndx_all;
            x5 = sp01_all;
            X = [x1, x2, x3, x4, x5];
           
        % ZNDX  
            Y = zndx_all;
            x1 = pdsi_all;
            x2 = phdi_all;
            x3 = pmdi_all;
            x4 = acres;
            x5 = sp01_all;
            X = [x1, x2, x3, x4, x5];
            
        
        % SP01  
            Y = sp01_all;
            x1 = pdsi_all;
            x2 = phdi_all;
            x3 = pmdi_all;
            x4 = acres;
            x5 = zndx_all;
            X = [x1, x2, x3, x4, x5];
            
        
    MulReg = regress(Y,X)
       
    scatter3(x1,x2,Y,'filled')
    hold on
    x1fit = min(x1):100:max(x1);
    x2fit = min(x2):10:max(x2);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    YFIT = MulReg(1) + MulReg(2)*X1FIT + MulReg(3)*X2FIT + MulReg(4)*X1FIT.*X2FIT;
    mesh(X1FIT,X2FIT,YFIT)
    xlabel('Weight')
    ylabel('Horsepower')
    zlabel('MPG')
    view(50,10)
    hold off
    
    
%}


%% T-Test
%     [h,pvalue,ci] = ttest(phdiPosfireNoZ(:,3),phdiPosfireNoZ(:,4),'Alpha',0.05)
%     [h,pvalue,ci] = ttest(acres,pdsi_all,'Alpha', 0.05)
%     [h,pvalue,ci] = ttest(acres,phdi_all,'Alpha',0.05)
%     [h,pvalue,ci] = ttest(phdi_pos(:,3),phdi_pos(:,4),'Alpha',0.05)
%     [h,pvalue,ci] = ttest(phdi_neg(:,3),phdi_neg(:,4),'Alpha',0.05)
%     [h,pvalue,ci] = ttest(acres,pmdi_all,'Alpha',0.05)
%     [h,pvalue,ci] = ttest(acres,zndx_all,'Alpha',0.05)
%     [h,pvalue,ci] = ttest(acres,sp01_all,'Alpha', 0.05)


%% Chi-Square
%    [h,pvalue,stats] = chi2gof(phdiPosfireNoZ(:,3),phdiPosfireNoZ(:,4),'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(acres,pdsi_all,'Alpha', 0.05)chi2gof
%     [h,pvalue,stats] = chi2gof(acres,phdi_all,'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(phdi_pos(:,3),phdi_pos(:,4),'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(phdi_neg(:,3),phdi_neg(:,4),'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(acres,pmdi_all,'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(acres,zndx_all,'Alpha',0.05)
%     [h,pvalue,stats] = chi2gof(acres,sp01_all,'Alpha', 0.05)


%% Sig. Testing
% t = tiledlayout(2,3);
% % figure()
% nexttile
% histogram(pdsi_all)
% %title('PDSI distribution (1984-2019)')
% title('PDSI (1984-2019)')
% xlabel('Index Value')
% ylabel('Count')
% %figure()
% nexttile
% histogram(phdi_all)
% %title('PHDI distribution (1984-2019)')
% title('PHDI (1984-2019)')
% xlabel('Index Value')
% ylabel('Count')
% %figure()
% nexttile
% histogram(pmdi_all)
% %title('PMDI distribution (1984-2019)')
% title('PMDI (1984-2019)')
% xlabel('Index Value')
% ylabel('Count')
% %figure()
% nexttile
% histogram(zndx_all)
% %title('ZNDX distribution (1984-2019)')
% title('ZNDX (1984-2019)')
% xlabel('Index Value')
% ylabel('Count')
% %figure()
% nexttile
% histogram(sp01_all)
% %title('SP01 distribution (1984-2019)')
% title('SP01 (1984-2019)')
% xlabel('Index Value')
% ylabel('Count')
% %figure()
% nexttile
% histogram(acres)
% %title('Burn Severity (acres) distribution (1984-2019)')
% title('Burn Severity (Acres)')
% ylabel('Acres'); set(gca, 'YScale', 'log')
% figure()
% histogram(fire_NoZ)
% title('Burn Severity (Acres) Distribution (Zero Values Removed) (1984-2019)')
% ylabel('Acres'); set(gca, 'YScale', 'log')