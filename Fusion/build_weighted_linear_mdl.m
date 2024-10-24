function [TIF_par,cluster] = build_weighted_linear_mdl(X,Y,band_codes_L,method,varargin)
%     close all;
% ks 2024/03/20: replace length() with size(X,1).
    p = inputParser;
    addParameter(p,'doplot',false);     % default if not plot obs and linear regression
    addParameter(p,'ir',[]);   
    addParameter(p,'ic',[]);       
    addParameter(p,'Band_plot',4); 
    addParameter(p,'cluster',[]);  % default is the first cluster
    addParameter(p,'d',[]);
    addParameter(p,'wfun',[]);
    
    % request user's input
    parse(p,varargin{:});
    doplot = p.Results.doplot;
    ir = p.Results.ir;
    ic = p.Results.ic;
    Band_plot = p.Results.Band_plot;
    cluster = p.Results.cluster;
    d = p.Results.d;
    wfun = p.Results.wfun;

    color_code = {'r','b','c','m'};
    

    %% Define weight functions
    switch wfun
        case 'Fair'
            w = @(x) 1./(1+abs(x));
            % tune = 1.4;
        case 'Cauchy'
            w = @(x) 1./(1+(x).^2);
            % tune = 2.385;
        case 'Sqrt'
            w = @(x) 1./(1+sqrt(x));
            % tune = 1;
    end

    %% Create lists to hold slope and intercept
    intercept = zeros(1,length(band_codes_L));
    if strcmp(method,'multi-variable')
        slope = zeros([length(band_codes_L) length(band_codes_L)]);
    else
        slope = zeros(1,length(band_codes_L));
    end
    Rsquared = zeros(1,length(band_codes_L));
    

    %% obtain the optimal k
    if ~isempty(cluster)
        k = max(unique(cluster(~isnan(cluster))));
    else
        % fprintf('Cluster info not given. Run TIF with all pairs.\r');
        k = 1;
        cluster = ones([size(X,1),1]);
    end
    % if more than one cluster, then find the centriods of each cluster
    if k>1   
        [cluster,C] = kmeans([X(:,:),Y(:,:)],k);  % update cluster!!
    else
        C = median([X(:,:),Y(:,:)],'omitnan'); % otherwise, use the median as centriods
    end
    
    %% Loop by k clusters
    for ik = 1:k

    TIF_par(ik).row = ir;
    TIF_par(ik).col = ic;
    TIF_par(ik).NumofObs = sum(~isnan(X(cluster==ik,1)));
    TIF_par(ik).cluster = cluster;


    if TIF_par(ik).NumofObs < 4 
        TIF_par(ik).QA = 0;
    else
        TIF_par(ik).QA = 1;
    end
%     fprintf('%d observation pairs. Run TIF. \r',TIF_par(ik).NumofObs);

    %% Loop run TIF by iband
    try
        for iband = band_codes_L 
            x = X(cluster==ik,iband);
            y = Y(cluster==ik,iband);
            r = d(cluster==ik);
            % r = (r-min(r))./(max(r)-min(r));   % standarize the distance (default is false)
            weight = w(r);
            % stem(weight)
            switch method
                case 'linear'  %% 1. Linear regression method/Ordinary least squares
                    mdl = fitlm(x,y,'Weights',weight);   % Add weighting func based on t-offset
                    % mdl = fitlm(x,y);   #No weighting func
                    intercept(iband) = mdl.Coefficients.Estimate(1);
                    slope(iband) = mdl.Coefficients.Estimate(2);
                    Rsquared(iband) = mdl.Rsquared.Ordinary;
                case 'robustfit'  %% 2. robust linear regression
                    % default robustfit weight function
                    % ks: use 'RobustOpts' to using the 'bisquare' weight function with the default tuning constant.
                    % mdl = fitlm(x,y,'RobustOpts',on); 
                    mdl = fitlm(x,y,'RobustOpts','on','Weights',weight); 
                    intercept(iband) = mdl.Coefficients.Estimate(1);
                    slope(iband) = mdl.Coefficients.Estimate(2);
                    Rsquared(iband) = mdl.Rsquared.Ordinary;
                case 'multi-variable-linear'
                    % mdl = fitlm(X(cluster==ik,:),y);
                    mdl = fitlm(X(cluster==ik,:),y,'Weights',weight);
                    intercept(iband) = mdl.Coefficients.Estimate(1);
                    slope(iband,:) = [mdl.Coefficients.Estimate(2:7)];   % row means response band, col means predictor bands 
                    Rsquared(iband) = mdl.Rsquared.Ordinary;
                case 'multi-variable-robustfit'
                    mdl = fitlm(X(cluster==ik,:),y,'RobustOpts','on','Weights',weight);
                    intercept(iband) = mdl.Coefficients.Estimate(1);
                    slope(iband,:) = [mdl.Coefficients.Estimate(2:7)];   % row means response band, col means predictor bands 
                    Rsquared(iband) = mdl.Rsquared.Ordinary;
            end  % end of switch method    
        end % end of iband    
        TIF_par(ik).Slopes = slope;
        TIF_par(ik).Intercepts = intercept;
        TIF_par(ik).Rsquared = Rsquared;
        TIF_par(ik).Centroid = reshape(C(ik,:),[6,2]);
    catch   % force TIF parameters to zeros if has errors, such as not enough obs pairs for robust fitting
        TIF_par(ik).QA = 0;
        if strcmp(method,'multi-variable')
            TIF_par(ik).Slopes = zeros(6,6);
        else
            TIF_par(ik).Slopes = zeros(1,6);
        end
        TIF_par(ik).Intercepts = zeros(1,6);
        TIF_par(ik).Rsquared = zeros(1,6);
        TIF_par(ik).Centroid = reshape(C(ik,:),[6,2]);
    end    % end of try
    % TIF_par

    %% plot matched obs and the linear regression (optional)
    if doplot 
        x = X(cluster==ik,Band_plot);
        y = Y(cluster==ik,Band_plot);
        ax1 = nexttile(4,[2,2]);
        switch ik
            case 1
                h0_1 = scatter(ax1,x,y,35,"o",'filled','MarkerEdgeColor','#c78af0','MarkerFaceColor','#c78af0');
                h0_1.DisplayName = 'Matching obs.Cluster 1';                   
            case 2
                h0_2 = scatter(ax1,x,y,35,"o",'filled','MarkerEdgeColor','#f1c232','MarkerFaceColor','#f1c232');
                h0_2.DisplayName = 'Matching obs.Cluster 2';
        end
        hold(ax1,'on')

        if Band_plot >=4  
            if max(y)>5000 && max(y)<8000
                xlim([0 6200])
                ylim([0 6200])
                xticks([0,1000,2000,3000,4000,5000,6000]);
                yticks([0,1000,2000,3000,4000,5000,6000]);
            
            elseif  max(y)>8000
                xlim([0 8200])
                ylim([0 8200])
                xticks([0,1000,2000,3000,4000,5000,6000,7000,8000]);
                yticks([0,1000,2000,3000,4000,5000,6000,7000,8000]);
            else
                xlim([0 4200])
                ylim([0 4200])
                xticks([0,1000,2000,3000,4000]);
                yticks([0,1000,2000,3000,4000]);
            end
          
        else
            xlim([0 3000])
            ylim([0 3000])
        end

%       h1: 1:1 line
        h1 = refline(1,0);  
        h1.Color = 'k';
        h1.DisplayName = '1:1 line';
        h1.LineWidth = 1.5;
        h1.LineStyle = "--";
        % h2: linear model
%         if TIF_par(ik).QA
%             lm = fitlm(x,y);
%             h2 = refline(lm.Coefficients.Estimate(2),lm.Coefficients.Estimate(1));
%         else
%             h2 = refline(0,0);
%         end
%         h2.Color = 	"#EDB120";
%         h2.DisplayName = 'OLS';
%         h2.LineWidth = 2;
%         h3: robustfit model 
        if TIF_par(ik).QA
            mdl_rf = fitlm(x,y,'RobustOpts','on','Weights',weight);   % weighted robustfit model
            % mdl_rf
            switch method
                case 'robustfit'      
                    h3 = refline(mdl_rf.Coefficients.Estimate(2),mdl_rf.Coefficients.Estimate(1));
                    if k==1    
                        h3.Color = 'k';
                    else
                        % h3.Color = 'g';
                        h3.Color = 	"#77AC30";
                    end
%                 otherwise
%                     continue;
                    % h3 = refline(mdl.Coefficients.Estimate(Band_plot+1),mdl.Coefficients.Estimate(1));
            end
            % % add R2 to the plot
            % R2 = TIF_par(ik).Rsquared;
            % if k==1
            %     text(4700,1500+500*(ik),['R$^{2}$ = ',num2str(R2(Band_plot))],'Interpreter', 'latex', 'FontSize', 12, 'Color','k');
            % else
            %     if ik==1
            %         text(4700,1500-500*(ik),['R$^{2}$ = ',num2str(R2(Band_plot))],'Interpreter', 'latex', 'FontSize', 12, 'Color','m');
            %     else
            %         text(4700,1500-500*(ik),['R$^{2}$ = ',num2str(R2(Band_plot))],'Interpreter', 'latex', 'FontSize', 12, 'Color','c');
            %     end
            % end
        else
            h3 = refline(0,0);
        end
        h3.DisplayName = 'robustfit';
        % h3.Color = 	"#77AC30";
        h3.LineWidth = 2;
%       h4: Centroid of clusters
        if k==1
            h4 = plot(C(ik,Band_plot),C(ik,Band_plot+6),'x','MarkerSize',15,'LineWidth',3,'Color','k'); 
            h4.DisplayName = 'Centroids.Cluster 1';
        else
            switch ik
                case 1
                    h4_1 = plot(C(ik,Band_plot),C(ik,Band_plot+6),'x','MarkerSize',15,'LineWidth',3,'Color',color_code{ik}); 
                    h4_1.DisplayName = 'Centroids.Cluster 1';
                case 2
                    h4_2 = plot(C(ik,Band_plot),C(ik,Band_plot+6),'x','MarkerSize',15,'LineWidth',3,'Color',color_code{ik}); 
                    h4_2.DisplayName = 'Centroids.Cluester 2';
            end

        end
        

        if Band_plot >=4
            if max(y)>5000 && max(y)<8000
                xlim([0 6200])
                ylim([0 6200])
                xticks([0,1000,2000,3000,4000,5000,6000]);
                yticks([0,1000,2000,3000,4000,5000,6000]);
            
            elseif  max(y)>8000
                xlim([0 8200])
                ylim([0 8200])
                xticks([0,1000,2000,3000,4000,5000,6000,7000,8000]);
                yticks([0,1000,2000,3000,4000,5000,6000,7000,8000]);
            else
                xlim([0 4200])
                ylim([0 4200])
                xticks([0,1000,2000,3000,4000]);
                yticks([0,1000,2000,3000,4000]);
            end
            
        else
            xlim([0 3000])
            ylim([0 3000])
        end


        xlabel('L30','FontSize',14)
        ylabel('S10','FontSize',14)
        title(['row/col=',num2str(ir),'/',num2str(ic),', Band ',num2str(Band_plot)]);
        if ik==1 && k==1
            legend([h0_1,h4],'Location','northwest');
        elseif ik==2 && k==2 
            legend([h0_1,h0_2,h4_1,h4_2],'Location','northwest');
        end
        fontname(gcf,"Serif");
%         set(gcf, 'color', 'none');   
    end % end of dolot
    end % end of ik   
    
end    % end of build_linear_mdl 

