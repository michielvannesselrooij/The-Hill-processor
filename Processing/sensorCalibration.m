function Fcal = sensorCalibration(F, showPlot)
% ------------------------------------------------------------------------
% Calibrate the input force with the calibration data
% INPUT:
% F         : Force data (vector)
% showPlot  : Show calibration data (1 or 0)
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    % Check input
    if ~exist('showPlot','var')
        showPlot = 0;
    end

    % Calibration data
%     refCal = [0,...
%               9.8100E-2,...
%               1.9600E-1,...
%               4.9000E-1,...
%               9.8100E-1,...
%               1.9610E+0];
%         
%     outCal = [-4.7732E-2,...
%               3.2641E-2,...
%               1.2333E-1,...
%               4.0021E-1,...
%               8.6602E-1,...
%               1.7996E+0];

    % New sensor [July, 2020]
    % No re-calibrations, used mV calibration from manufacturer in LabView
    refCal = [0, 2];
    outCal = [0, 2];
        
    % Calibrate
    [p, ~, mu] = polyfit(outCal, refCal, 1);
    Fcal       = polyval(p, F, [], mu); 
        
    % Show calibration (optional)
    if showPlot
       figure; hold on; box on; grid minor;
       title('Force data calibration');
       xlabel('Uncalibrated [N]');
       ylabel('Reference [N]');
       
       xq = linspace(min(outCal)*1.2, max(outCal)*1.2, 200);
       yq = polyval(p, xq, [], mu);
       plot(xq,yq,'k-');
       
       plot(outCal, refCal, 'ks');  % Calibration data
       plot(F, Fcal,'ro');          % Measurement data
       
    end

end