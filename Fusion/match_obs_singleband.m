function [X,Y] = match_obs_singleband(clrx_L30,clrx_S2,clry_L30,clry_S2,band_codes_L,band_codes_S,t_threshold,hide_date,match_method)
% find matched observations and return X,Y to build linear regression.
%   X - matched L30 SR in all bands         (predictor variable)
%   Y - matched Sentinel-2 SR in six bands  (response variable)
% Inputs:
% t_threshold - The temporal threshold used to match L30 and S2 observations.
%            ex: +-1 day, +-2 days
%
    switch match_method
        case 'individual'
            X = NaN([3000,1]);  % L30 surface reflectance, the size(1) of X should be large enough. Here I give 3000.
            Y = NaN([3000,1]);  % S10 surface reflectance 
            k = 1;
        otherwise
            match_id = zeros([length(clrx_L30),4]);  % matched observation ids
            X = NaN([length(clrx_L30),1]);  % L30 surface reflectance
            Y = NaN([length(clrx_L30),1]);  % S10 surface reflectance 
    end

    for i = 1:length(clrx_L30)  
        idate = clrx_L30(i);
        % ex. To test resuls, we hide 20210616 for reference
        if idate == datenum(hide_date)  
            continue;
        end
        t_dif = abs(idate-clrx_S2); % temporal difference (day)
        if min(t_dif)<t_threshold
            close_id = find(t_dif==min(t_dif));  % The best: L30 and S10 at the same day!
        else
            close_id = find(t_dif<=t_threshold); % Sacrifize temporal accuracy for matchness.
        end
        if isempty(close_id)   % no matched observation, give nan
            close_id = nan;
        else                    % matched obs exist
            for iband = band_codes_L
                switch match_method
                    case 'individual'
                        % new code (20230321): allow one L30 to match
                        % multiple S10 - means more pairs for each L30.
                        Y(k:k+length(close_id)-1,iband) = clry_S2(close_id,band_codes_S(iband));
                        X(k:k+length(close_id)-1,iband) = clry_L30(i,iband); % 2: Green, 3:Red, 4:NIR, 5: SWIR1, 6: SWIR2 
                    case 'average'
                        %old code: one L30 ONLY match one S10 (the average value)
                        if length(close_id)>1   % more than one matched obs and none at the same day (eg. 1 2, 2 2), give the mean value
                            match_id(i,1:length(close_id))= close_id;
                            sr_S2 = mean(clry_S2(close_id,band_codes_S(iband)));  % 2: Green, 10:NNIR
                        else
                            match_id(i,1) = close_id;
                            sr_S2 = clry_S2(close_id,band_codes_S(iband));
                        end  
                        Y(i,iband) = sr_S2;
                        X(i,iband) = clry_L30(i,iband); % 2: Green, 3:Red, 4:NIR, 5: SWIR1, 6: SWIR2
                    case 'first'
                        % use the earlist S10 obs to match L30 (the first value).
                        if length(close_id)>1   % more than one matched obs and none at the same day (eg. 1 2, 2 2), give the mean value
                            close_id = min(close_id);
                            sr_S2 = clry_S2(close_id,band_codes_S);  % 2: Green, 10:NNIR
                        else
                            match_id(i,1) = close_id;
                            sr_S2 = clry_S2(close_id,band_codes_S);
                        end  
                        Y(i) = sr_S2;
                        X(i) = clry_L30(i); % 2: Green, 3:Red, 4:NIR, 5: SWIR1, 6: SWIR2
                end   % end of switch
            end % end of iband
            if strcmp(match_method,'individual')
                k = k+length(close_id);
            end
        end % end of isempty(close_id)
    end % end of length(clrx_L30)
end


