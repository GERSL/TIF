function plotTimeSeries(clrx_S10,clry_S10,clrx_L30,clry_L30,Band_plot,pt_lat,pt_lon,varargin)
% plot original clear L30 and S10 time series.
% plot TIF predictions.
%
% ks 20230817: add line 58-69 to use centroid to determine which TIF par
% for prediction the hide_L30.
    addpath(fullfile(pwd, 'Fusion'));
    p = inputParser;
    addParameter(p,'plot_pred',false);
    addParameter(p,'TIF_coefficient',[]);
    addParameter(p,'test_dates',[]);
    addParameter(p,'cluster',[]);
    addParameter(p,'multivariable',false);

    parse(p,varargin{:});
    plot_pred = p.Results.plot_pred;
    TIF_coefficient = p.Results.TIF_coefficient;
    test_dates = p.Results.test_dates;
    cluster = p.Results.cluster;
    multivariable = p.Results.multivariable;


    % set(gcf,'Position',[100 100 800 800]);
    set(gcf,'Position',[50 100 1100 400]);
    set(gca,'FontSize',13);
    % set(gcf, 'Color', 'none');
    set(gcf, 'Color', 'w');

    if ~plot_pred
        tiledlayout(2,5,"TileSpacing","compact","Padding","compact");     
    end
    ax0 = nexttile(6,[1,3]);

    % band codes for Landsat and Sentinel-2 
    band_codes_L = [1,2,3,4,5,6];
    band_codes_S = [1,2,3,4,569];
    daterange =[datenum(2013,1,1), datenum(2022,12,31)];


    if isempty(cluster)
        cluster = ones([length(clrx_L30),1]);
    end
    
    % find the first nonNaN idx
    nnan = find(~isnan(cluster));
    first_nnan = nnan(1);
    new_cluster = NaN(size(cluster));
    % fill NaN in idx
    for i = 1:length(cluster)
        if isnan(cluster(i))
            value = clry_L30(i,1:6);
            k = length(TIF_coefficient);
            if k>1 % if more than one TIF outcomes, determine which TIF to use based on the TIF_coefficient.Centroid
               for j = 1:k
                   tmp = TIF_coefficient(j).Centroid;
                   point = tmp(:,1)';
                   d(j) = pdist([value;point],'euclidean');
               end
               new_cluster(i) = find(d==min(d));
            else
               new_cluster(i) = 1;
            end   % end of if ik>1
        else
            new_cluster(i) = cluster(i);
        end   % end of if isnan(cluster(i))
    end    % end of for i = 1:length(cluster)

    % Plot raw L30 and S10 time series
    figure(1)
    % doy = clrx_S10 - datenum(2013,1,1);
    p1 = plot(clrx_S10,clry_S10(:,band_codes_S(Band_plot)), 'o','MarkerEdgeColor','#f78b8b','MarkerFaceColor','#f78b8b', 'Markersize', 4,'DisplayName', 'S10');
    hold(ax0,'on')
    % doy = clrx_L30 - datenum(2013,1,1);
    p2 = plot(clrx_L30,clry_L30(:,band_codes_L(Band_plot)), 'o', 'MarkerEdgeColor','#8bb8f7','MarkerFaceColor','#8bb8f7', 'Markersize', 4,'DisplayName', 'L30');
    lgd = legend([p1,p2],'Location','best');
    
    % Loop all TIF relationships and plot predictions
    if plot_pred
        % create two lists to save accuracy merics
        rmse = zeros(length(TIF_coefficient),1);
        AD = zeros(length(TIF_coefficient),1);
        reference = [];
        prediction = [];
        for ik = 1:length(TIF_coefficient)
            if TIF_coefficient(ik).QA 
                if multivariable
                    try
                        slope_iband = TIF_coefficient(ik).Slopes(Band_plot,:);
                    catch
                        slope_iband = TIF_coefficient(ik).Slopes(1,:);
                    end
                    a = clry_L30(new_cluster==ik,:);
                    b = slope_iband';
                    pred = a*b+TIF_coefficient(ik).Intercepts(Band_plot);
                else
                    pred = clry_L30(new_cluster==ik,:).*TIF_coefficient(ik).Slopes+TIF_coefficient(ik).Intercepts;
                end
                % Adjust values exceed [0,10000]
                pred(pred<0) = 0;
                pred(pred>10000) = 10000;
                % Calculate RMSE and AD
                t_threshold = 1;
                hide_date = [];
                % match single babd prediction and reference
                if multivariable
                    [pred_sb,ref_sb] = match_obs_singleband(clrx_L30(new_cluster==ik),clrx_S10,pred,clry_S10,band_codes_L(Band_plot),band_codes_S(Band_plot),t_threshold,hide_date,'first');
                else
                    [pred_sb,ref_sb] = match_obs_singleband(clrx_L30(new_cluster==ik),clrx_S10,pred(:,Band_plot),clry_S10,band_codes_L(Band_plot),band_codes_S(Band_plot),t_threshold,hide_date,'first');
                end
                reference = [reference;ref_sb];
                prediction = [prediction;pred_sb];
                
                % display TIF predictions on the time series plot
                if length(TIF_coefficient)>1
                    switch ik
                        case 1
                            if multivariable
                                p3_1 = plot(clrx_L30(new_cluster==ik),pred(:), 'o', 'Color','#c78af0', 'LineWidth',1, 'Markersize', 8,'DisplayName', ['TIFpred.Cluster ',int2str(ik)]);
                            else
                                p3_1 = plot(clrx_L30(new_cluster==ik),pred(:,band_codes_L(Band_plot)), 'o', 'Color','#c78af0','LineWidth',0.8, 'Markersize', 8,'DisplayName', ['TIFpred.Cluster ',int2str(ik)]);
                            end
                        case 2
                            if multivariable
                                p3_2 = plot(clrx_L30(new_cluster==ik),pred(:), 'o', 'Color','#f1c232','LineWidth',1, 'Markersize', 8,'DisplayName', ['TIFpred.Cluster ',int2str(ik)]);
                            else
                                p3_2 = plot(clrx_L30(new_cluster==ik),pred(:,band_codes_L(Band_plot)), 'o', 'Color','#f1c232','LineWidth',0.8, 'Markersize', 8,'DisplayName', ['TIFpred.Cluster ',int2str(ik)]);
                            end
                    end   % end of switch ik
                else
                    p3_1 = plot(clrx_L30(new_cluster==1),pred(:,band_codes_L(Band_plot)), 'o', 'Color','k', 'LineWidth',0.3,'Markersize', 8,'DisplayName', 'TIFpred.Cluster 1');
                end
            end  % end of TIF_coefficient(ik).QA

        end   % end of ik = 1:length(TIF_coefficient)
        % rmse = CalRMSE(reference(~isnan(reference)),prediction(~isnan(reference)));
        % AD = CalBias(reference(~isnan(reference)),prediction(~isnan(reference)));
        % fprintf('RMSE = %.3f\n',rmse);
        % fprintf('AD = %.3f\n',AD);

        if length(TIF_coefficient)==1
            if TIF_coefficient(1).QA
                % lgd = legend([p1,p2,p3_1],'Location','bestoutside');
                lgd = legend([p1,p2,p3_1],'Location','northwest');
            end

        else
            if TIF_coefficient(1).QA && TIF_coefficient(2).QA
                % lgd = legend([p1,p2,p3_1,p3_2],'Location','bestoutside');
                lgd = legend([p1,p2,p3_1,p3_2],'Location','northwest');
            end
        end  % end of if length(TIF_coefficient)==1
    end   % end of if plot_pred

    % Plot the test dates as dotted lines.
    if ~isempty(test_dates)
        xline(datenum(test_dates),'--k','HandleVisibility','off');
%         xline(datenum(test_dates{length(test_dates)}),'--k','DisplayName','test dates');
    end
    
    
    % datetick('x', 2, 'keeplimits');
    % datetick('x', 10, 'keeplimits');

    ax = gca;
    % set x-axis range
    xlim([datenum('2013-01-01'),datenum('2022-12-31')]);
    datetick('x', 10, 'keeplimits');
    ax.FontSize = 12;
    if Band_plot==4
        % ylim([0 5500]);
        ylim([0,8000]);
        % ylim([1000 4500]);
        % ylim([1000, 6500])
    else
        ylim([0,6000]);
        % ylim([1000,6000]);
    end
    ylabel(['Band ',num2str(Band_plot)],'FontSize',16);
    if pt_lat~=0
        title(sprintf('Lat/Lon: %.4f %.4f', pt_lat(1), pt_lon(1)),'FontSize',16);
    end

    fontname(gcf,"Lucida Bright")

    %% Plot the temporal distribution of clusters (optional)
%     figure()
%     set(gcf,'Position',[100 100 800 300]);
%     stem(clrx_L30,new_cluster,':r','DisplayName','all obs')
%     hold on
%     stem(clrx_L30,cluster,':b','DisplayName','matching obs')
%     hold on
%     if ~isempty(test_dates)
%         for i =1:length(test_dates)
%             stem(datenum(test_dates),1,'diamondk','filled','DisplayName','test date');
% %             stem(clrx_L30(clrx_L30==datenum(test_dates)),new_cluster(clrx_L30==datenum(test_dates)),'diamondk','filled','DisplayName','test date')
% %             stem(clrx_L30(clrx_L30==datenum(test_dates{i})),new_cluster(clrx_L30==datenum(test_dates{i})),'diamondk','filled','DisplayName','test date')
%         end
%     end
%     ylim([0 2.5])
%     yticks([0,1,2])
%     xlim([datenum('2013-01-01'),datenum('2021-12-31')]);
%     datetick('x', 10, 'keeplimits');
%     legend('Location','bestoutside')
% 
%     fontname(gcf,"Lucida Bright")
% 


end


function RMSE = CalRMSE(Ref, Pred)
    dif(:) = Ref(:) - Pred(:);
    dif(:) = dif(:).^2;
    RMSE = sqrt(mean2(dif(:)));      
end

function [AAD, AD] = CalBias(Ref, Pred)

    AAD = mean(abs(Pred-Ref));
    AD = mean(Pred-Ref);

end