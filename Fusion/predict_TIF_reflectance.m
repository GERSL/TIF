function prediction = predict_TIF_reflectance(clry_L, TIF_par, bands, multi_variable)
    value = clry_L(1:6);
    k = length(TIF_par);
    % determine TIF cluster
    if k>1 % if more than one TIF outcomes, determine which TIF to use based on the TIF_par.Centroid
       for j = 1:k
           tmp = TIF_par(j).Centroid;
           point = tmp(:,1)';
           d(j) = pdist([value;point],'euclidean');
       end
       cluster = find(d==min(d),1);  % return only one value when there are same d
    else
       cluster = 1;
    end   % end of if ik>1
    
    % calculate TIFprediction
    prediction = NaN(length(bands),1);
    if multi_variable
        for band_id = 1: length(bands)
            band = bands(band_id);
            if TIF_par(cluster).QA 
                try
                    slope_iband = TIF_par(cluster).Slopes(band,:);
                catch
                    slope_iband = TIF_par(cluster).Slopes(1,:);
                end
                a = clry_L;
                b = slope_iband';
                pred = a*b+TIF_par(cluster).Intercepts(band);
                prediction(band_id) = pred;
            end
        end   % end of band_id
    else
        if TIF_par(cluster).QA 
            prediction = clry_L.*TIF_par(cluster).Slopes+TIF_par(cluster).Intercepts;
        end
    end  % end of if multi_variable
end   % end of func