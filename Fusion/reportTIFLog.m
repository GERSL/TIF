function reportTIFLog(folderout, ntasks, regress_method, wfun,  t_threshold)
%REPORTLOG This is to report the log of TIF.
 % record mutilple logs ...
	TIF_v = 1.00; % version 1.00
    HLS_v = 2.0;   % version v2.0
    % fileID = fopen(fullfile(folderout, 'TIF_log.txt'),'a'); % Open or create new file for writing. Append data to the end of the file.
    fileID = fopen(fullfile(folderout, 'TIF_log.txt'),'w');   
    % fileID = fopen('COLD_log.txt','w'); % Open or create new file for writing. Discard existing contents, if any.

    % record the starting time of processing.
    % fprintf(fileID,'Hide date = %s\r\n', hidedate); 
    fprintf(fileID,'Starting Time = %s\r\n', datetime('now')); 
    fprintf(fileID,'Total Cores = %d\r\n', ntasks); 

    % % analysis scale
    % fprintf(fileID,'Analysis scale = %s\r\n', analysis_scale);

    % regression method
    fprintf(fileID,'Regression method = %s\r\n', regress_method);
    % observation matching temporal window
    fprintf(fileID,'Observation matching window = %d day\r\n',t_threshold);
    % observation matching temporal window
    fprintf(fileID,'Observation weight function = %s \r\n',wfun);


    % TIF Version
    fprintf(fileID,'TIF Version = %.2f\r\n',TIF_v);
    % HLS Version
    fprintf(fileID,'HLS Data Version = %.2f\r\n',HLS_v);

    % updates
    fprintf(fileID,'******************************************************************************************************\r\n');
    fprintf(fileID,'Revisions: $ Date: 10/24/2024 $ Copyright: GERS Lab, UCONN\r\n');
    
    fprintf(fileID,'Version 1.00   Time-series-based Image Fusion (10/24/2024) \r\n');
    fprintf(fileID,'******************************************************************************************************\r\n');

    fclose(fileID);
end

