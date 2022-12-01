%=========================================================================
% Description: The code below performs time series analysis on Palmer and Standard Precipitiation Drought
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
%       FIGURES - Output time series
% 
% Code: 
%%      This is based on HW2_Solutions.m by Dr. Ruixin Yang (George Mason University)
%           -originally sourced from code by John H. Powell (PlotNino.m) in 2009
%           and modified by Dr. Ruixin Yang (D:\Yang\Work-Research\Code\examples\read_ASCII
%           folder).
%   
% Author: Brendan Casey
% Assign: Graduate Thesis (M.S. in Geographic and Cartographic Sciences)
% Defense Date: 20220818
% Modification Records: 20221130 (cleaned code, added descrip.)
%                       20221201 (cleaned code, rev. descrip.)
%=========================================================================


%% pick the data for given year period only.
start_year = 1984;
start_year_off = 1983;
end_year = 2019;
end_year_off = 2018;
n_month = 12*(end_year - start_year +1);
n_month_off = 12*(end_year_off - start_year_off +1);
wfire = fire(1:n_month, :);

               
%% plot whole time series in one plot. 

% Burn
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('12-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, acres);
    title_str =sprintf('Burn Severity for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time (mm/yy)');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Burn Severity (Total Acres)');
    grid on;

% PDSI
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('11-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, pdsi_all);
    title_str =sprintf('PDSI for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Palmer Index Values');
    grid on;
% PHDI
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('12-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, phdi_all);
    title_str =sprintf('PHDI for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time (mm/yy)');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Palmer Hydro Index Values');
    grid on;

% PMDI
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('12-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, pmdi_all);
    title_str =sprintf('PMDI for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time (mm/yy)');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Palmer Modified Index Values');
    grid on;

% ZNDX
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('12-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, zndx_all);
    title_str =sprintf('ZNDX for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time (mm/yy)');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Palmer-Z Index Values');
    grid on;

% SP01
figure; 
    start_str = sprintf('01-01-%d', start_year);
    end_str = sprintf('12-31-%d', end_year);
    startDate = datenum(start_str); 
    endDate = datenum(end_str);
    xData = linspace(startDate,endDate,n_month);
    plot(xData, sp01_all);
    title_str =sprintf('SP01 for Jan%d-Dec%d Period',start_year, end_year);
    title(title_str);
    xlabel('Time (mm/yy)');
    ax = gca;
    ax.XTick = xData;
    datetick('x','mm/yy');
    ylabel('Std. Prec. Index (one-month) Values');
    grid on;


% Plot the same time series with subplots to better illustrate 
% values 

% Burn
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 1; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),acres(Ind));
    ax = gca;
    ax.XTick = xData(Ind);
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('Burn Severity Time Series, Jan%d-Dec%d',start_year, end_year);
title(title_str);


% PDSI
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 6; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),pdsi_all(Ind));
    %2nd yaxis (wildfire)
        yyaxis right
        plot(xData(Ind), acres(Ind), 'k');
        xlabel('Year');
        ystring = sprintf('TFP for %s',acres);
        ylabel(ystring);
    ax = gca;
    ax.XTick = xData(Ind);
%     set(gca,'XTick',1:x_tick_interval:x_length)
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
    yyaxis left
    ylabel('Index Values');
    %legend('PDSI','Acres')
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('PDSI Time Series with Total Burn Acreage, Jan%d-Dec%d',start_year, end_year);
title(title_str)


% PHDI
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 6; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),phdi_all(Ind));
    %2nd yaxis (wildfire)
        yyaxis right
        plot(xData(Ind), acres(Ind), 'k');
        xlabel('Year');
        ystring = sprintf('TFP for %s',acres);
        ylabel(ystring);
    ax = gca;
    ax.XTick = xData(Ind);
%     set(gca,'XTick',1:x_tick_interval:x_length)
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
    yyaxis left
    ylabel('Index Values');
    %legend('PHDI','Acres')
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('PHDI Time Series with total burn acreage, Jan%d-Dec%d',start_year, end_year);
title(title_str)


% PMDI
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 6; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),pmdi_all(Ind));
    %2nd yaxis (wildfire)
        yyaxis right
        plot(xData(Ind), acres(Ind), 'k');
        xlabel('Year');
        ystring = sprintf('TFP for %s',acres);
        ylabel(ystring);
    ax = gca;
    ax.XTick = xData(Ind);
%     set(gca,'XTick',1:x_tick_interval:x_length)
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
    yyaxis left
    ylabel('Index Values');
    %legend('PMDI','Acres')
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('PMDI Time Series with Total Burn Acreage, Jan%d-Dec%d',start_year, end_year);
title(title_str)


% ZNDX
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 6; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),zndx_all(Ind));
    %2nd yaxis (wildfire)
        yyaxis right
        plot(xData(Ind), acres(Ind), 'k');
        xlabel('Year');
        ystring = sprintf('TFP for %s',acres);
        ylabel(ystring);
    ax = gca;
    ax.XTick = xData(Ind);
%     set(gca,'XTick',1:x_tick_interval:x_length)
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
    yyaxis left
    ylabel('Index Values');
    %legend('ZNDX','Acres')
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('ZNDX Time Series with Total Burn Acreage, Jan%d-Dec%d',start_year, end_year);
title(title_str)


% SP01
figure
n_subplot = 3; 
x_length = 144; % number of months in one subplot
x_tick_interval = 6; %Plot a tick every six months

for i=1:n_subplot
    subplot(n_subplot,1,i);
    istart = x_length*(i-1)+1;
    iend = x_length*i;
    if(iend>n_month) iend = n_month; end
    Ind=(istart:iend);
    plot(xData(Ind),sp01_all(Ind));
    %2nd yaxis (wildfire)
        yyaxis right
        plot(xData(Ind), acres(Ind), 'k');
        xlabel('Year');
        ystring = sprintf('TFP for %s',acres);
        ylabel(ystring);
    ax = gca;
    ax.XTick = xData(Ind);
%     set(gca,'XTick',1:x_tick_interval:x_length)
    datetick('x','mm/yy');
    ylabel('Acres');
    grid on;
    yyaxis left
    ylabel('Index Values');
    %legend('SP01','Acres')
end
xlabel('Date');
subplot(n_subplot,1,1);
title_str =sprintf('SP01 Time Series with Total Burn Acreage, Jan%d-Dec%d',start_year, end_year);
title(title_str)


% Compute monthly climo & anomalies
for i=1:12
    ind = find(fire(:,2)==i); % find data indices in the original time 
                                % series for a given month
    MonClimo(i)=mean(acres(ind));
    NinoAnomCalc(ind)=acres(ind)-MonClimo(i);
end


% Plot means
figure
xDatam = linspace(datenum('01-01-2009'),datenum('12-31-2009'),12);
plot(xDatam, MonClimo);
title('Burn Severity Monthly Climatology (1984-2019)');
xlabel('Months');
ax = gca;
ax.XTick = xDatam;
datetick('x','mmm','keepticks');
ylabel('Average Burn Severity (Acres)');
grid on;


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T-Test 
    % Burn Acres
    [h,p,ci,stats] = ttest(acres,0,'Alpha',0.05)
    % PDSI
    [h,p,ci,stats] = ttest(pdsi_all,0,'Alpha',0.05)
    % PHDI
    [h,p,ci,stats] = ttest(phdi_all,0,'Alpha',0.05)
    % PMDI
    [h,p,ci,stats] = ttest(pmdi_all,0,'Alpha',0.05)
    % ZNDX
    [h,p,ci,stats] = ttest(zndx_all,0,'Alpha',0.05)
    % SP01
    [h,p,ci,stats] = ttest(sp01_all,0,'Alpha',0.05)




% F-Test 
    % Burn Acres
    [h,p,ci,stats] = vartest(acres,0,'Alpha',0.05)
    % PDSI
    [h,p,ci,stats] = vartest(pdsi_all,0,'Alpha',0.05)
    % PHDI
    [h,p,ci,stats] = vartest(phdi_all,0,'Alpha',0.05)
    % PMDI
    [h,p,ci,stats] = vartest(pmdi_all,0,'Alpha',0.05)
    % ZNDX
    [h,p,ci,stats] = vartest(zndx_all,0,'Alpha',0.05)
    % SP01
    [h,p,ci,stats] = vartest(sp01_all,0,'Alpha',0.05)
%}
% Scatter plots
    t = tiledlayout(3,2);
    % Burn vs PDSI
    %figure;
    nexttile
    scatter(pdsi_all, acres,2);
    ylabel('Acres Burned');
    xlabel('PDSI Values');
    grid on;
    % Burn vs PHDI
    %figure;
    nexttile
    scatter(phdi_all, acres,2);
    ylabel('Acres Burned');
    xlabel('PHDI Values');
    grid on;
    % Burn vs PMDI
    %figure;
    nexttile
    scatter(pmdi_all, acres,2);
    ylabel('Acres Burned');
    xlabel('PMDI Values');
    grid on;
    % Burn vs ZNDX
    %figure;
    nexttile
    scatter(zndx_all, acres,2);
    ylabel('Acres Burned');
    xlabel('ZNDX Values');
    grid on;
    % Burn vs SP01
    %figure;
    nexttile
    scatter(sp01_all, acres,2);
    ylabel('Acres Burned');
    xlabel('SP01 Values');
    grid on;
