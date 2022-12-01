% Burn Severity:
    % Categories:
        % All:
            fire = load('wildfire_num_only.txt'); 
            sz_f = size(fire);
            rows_f = sz_f(1);
            cols_f = sz_f(2);
            %time_f = (1:1:((rows_f)*(cols_f-2)))';

     

        % Burn Levels:
            %Create variables (single column) for burn severity and years.

            % Severity (all years)
                unburn = fire(:,3);
                low = fire(:,4);
                mod = fire (:,5);
                hi = fire (:,6);
                inc_gre = fire(:,7);
                non_pr = fire(:,8);
                event = fire(:,9); %total fire events (days)
                acres = fire(:,80); %total acres burned (low+moderate+high...etc)
            

   
% Drought:
    pdsi = load('pdsi.txt'); 
    pdsi_offset = load('pdsi_offset.txt');
    phdi = load('phdi.txt');
    phdi_offset = load('phdi_offset.txt');
    pmdi = load('pmdi.txt');
    pmdi_offset = load('pmdi_offset.txt');
    zndx = load('zndx.txt');
    zndx_offset = load('zndx_offset.txt');
    sp01 = load('sp01.txt');
    sp01_offset = load('sp01_offset.txt');
    

    % Palmer Drought Severity Index
        sz_pdsi = size(pdsi);
        rows_pdsi = sz_pdsi(1);
        cols_pdsi = sz_pdsi(2);
        time_pdsi = (1:1:((rows_pdsi)*(cols_pdsi-2)))';

    % Palmer Hydrologic Drought Index
        sz_phdi = size(phdi);
        rows_phdi = sz_phdi(1);
        cols_phdi = sz_phdi(2);
        time_phdi = (1:1:((rows_phdi)*(cols_phdi-2)))';
   
    % Modified Palmer Drought Severity Index
        sz_pmdi = size(pmdi);
        rows_pmdi = sz_pmdi(1);
        cols_pmdi = sz_pmdi(2);
        time_pmdi = (1:1:((rows_pmdi)*(cols_pmdi-2)))';

    % Palmer "Z" Index
        sz_zndx = size(zndx);
        rows_zndx = sz_zndx(1);
        cols_zndx = sz_zndx(2);
        time_zndx = (1:1:((rows_zndx)*(cols_zndx-2)))';

    % Standardized Precipitation Index (one-month)
        sz_sp01 = size(sp01);
        rows_sp01 = sz_sp01(1);
        cols_sp01 = sz_sp01(2);
        time_sp01 = (1:1:((rows_sp01)*(cols_sp01-2)))';


%%Create variables (single column) for drought months
        
% Palmer Drought Severity Index
            pdsi_all = pdsi(:,2:13);
            pdsi_all = reshape(pdsi_all',1,numel(pdsi_all));
            pdsi_all = pdsi_all(:);
  
% Palmer Hydrologic Drought Index
            phdi_all = phdi(:,2:13);
            phdi_all = reshape(phdi_all',1,numel(phdi_all));
            phdi_all = phdi_all(:);

% Modified Palmer Drought Severity Index
            pmdi_all = pmdi(:,2:13);
            pmdi_all = reshape(pmdi_all',1,numel(pmdi_all));
            pmdi_all = pmdi_all(:);

% Palmer "Z" Index
            zndx_all = zndx(:,2:13);
            zndx_all = reshape(zndx_all',1,numel(zndx_all));
            zndx_all = zndx_all(:);
            
% Standardized Precipitation Index (one month)
            sp01_all = sp01(:,2:13);
            sp01_all = reshape(sp01_all',1,numel(sp01_all));
            sp01_all = sp01_all(:);
            

         