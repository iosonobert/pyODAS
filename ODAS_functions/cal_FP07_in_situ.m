%% cal_FP07_in_situ
% Calibrate thermistor probe using in situ data
%%
% <latex>\index{Functions!cal\_FP07\_in\_situ}</latex>
%
%%% Syntax
%   [T_0, beta, Lag] = cal_FP07_in_situ( file_name, T_ref_string, T_string, SN, ... )
%
% * [file_name] Name of mat-file containing data used to
%       calibrate a thermistor.
% * [T_ref_string] Name of vector within the mat-file that contains the
%       reference temperature in degrees celsius. Usually from a SBE4F
%       thermometer or a JAC_CT.
% * [T_string] Name of thermistor to calibrate, typically 'T1' or 'T2'.
% * [SN] Serial number of thermistor.
% * [...] Optional elements describing a profile.  Can be provided in a
%       structure or as key / value pairs.  See below for more details.
% * []
% * [T_0] Value of parameter T0, used in the Steinhart-Hart equation. When
%       called with no input parameters, a structure containing default 
%       input parameter values is returned for reference.
% * [beta] beta coefficients, in ascending order, of the fit to the SH 
%       equation. i.e. beta_1, beta_2, beta_3
% * [Lag] Delay in seconds between the thermistor and the reference
%       thermometer. Typically a negative value because the reference
%       sensor is usually behind the thermistor being calibrated.
%
%%% Description
%
% Function to calibrate a FP-07 thermistor probe using in-situ data.  Reliable
% temperature data must be within the specified mat-file. They will be the 
% refercence temperature -- usually a Sea-Bird SBE3F or a JAC-CT. 
%
% This function makes 5 figures. Figure 1 shows you the portion of the data
% file that will be used for the calibration. It then detrends the
% thermistor data and scales it so that it is approximately aligned with
% the reference thermometer (figure 2). It then plots the cross-correlation
% coefficient between the thermistor and the reference thermometer and
% estimates the lag between these two signals (Figure 3). Next it plots the
% the natural logarithm of the resistance ratio against the inverse of the
% absolute temperature to show the quality of the regression (Figure 4).
% Finally, it plots a depth profile of the lagged reference temperature,
% the newly calibrated thermistor temperature, and their difference.
%
% Optional input parameters include;
%
% * [profile_num] Integer identifying the profile number to be
%       analyzed. Default = 1. Some files may contain multiple profiles.
% * [vehicle_info] A structure found in all *.mat data files that is used
%       to determine the direction of profiling. If empty, the direction
%       defaults to 'down'. For glider data, the functions looks at up- and
%       down-profiles. Default = [].
%      'up' or 'down' - default = 'down'.
% * [profile_min_duration] Minimum time, in seconds, for a profile to be deemed a
%       profile.  Direction of travel must be monotonic, rate of change of 
%       pressure must be above a minimum level, and pressure must 
%       exceed a minimum value. Default = 20 s.
% * [profile_min_P] Minimum pressure for a profile. Default = 1 dBar
% * [profile_min_W] Minimum vertical speed of a profile, actually |dP/dt|.
%       Default = 0.4 dBar/s. Use a smaller value, ~0.1 dBar/s, for
%       gliders.
% * [order] Fit order to the Steinhart-Hart equation. Value can be 1,
%       2, or 3. Default = 2.

%___________________________
%
% Version History
%
% * 2013-12-05 (RGL) original version.
% * 2013-12-06 (RGL) added varargin to define a profile with default values.
% * 2015-04-10 (WID) revised documentation for publishing.
% * 2015-04-27 (RGL) modified to allow the specification of the fit order.
% * 2015-07-27 (WID) use texstr in place of fix_underscore.
% * 2015-07-29 (WID) return default values when called with no input
%                       parameters.
% * 2015-10-27 (RGL) Changed description section.
% * 2016-06-07 (RGL) Changed legend call, added clf to the start of figures.
% * 2016-11-10 (RGL) Changed the call to the deconvolve function so that it
%      includes the thermistor signal without pre-emphasis. odas_p2mat uses
%      the the thermistor without pre-emphasis (if it is available) for
%      conversion in to physical units. 
%      Consequently, the thermistor signal without pre-emphasis MUST be
%      called with this in situ calibration function. Otherwise, there is a
%      substantial error due to the offset (order ~10 counts) that may be
%      present in the signal with pre-emphasis. This version also includes
%      a test for the existence of the signal without pre-emphasis so that
%      it does not bomb in its absence. I also beautified the legends and
%      labels.
% * 2016-11-10 (RGL) Added a low-pass filter to the thermistor data, in the
%      case of a JAC-T reference, in order to get a tighter regression. The
%      cut-off frequency is sped dependent and follows the recommendation
%      for calculating salinity. There is an improvement when the
%      thermistor is filtered.
% * 2016-11-15 (RGL) Added fc in case of Sea-Bird thermometer.
% * 2017-11-28 (RGL) Added ability to handle both type=therm abd type=t_ms.
% * 2017-11-30 (RGL) Some more display changes and warning will be
%      suppressed in case the channel without pre-emphasis does not exist.
%      Made the naming of parameters (fields) consistent with usage in
%      quick_look.


function [T_0,beta,Lag] = cal_FP07_in_situ(file_name,T_ref_string,T_string,SN,varargin)

%%%%
% Default values for optional fields
default_vehicle_info         = []; % will trigger profile_dir = 'down'
default_profile_min_P        = 1; % in dBar
default_profile_min_W        = 0.4; % in dbar/s
default_profile_num          = 1; % process the first profile
default_profile_min_duration = 20; % minimum duration [s] for a segment to be considered a profile
default_order                = 2; % The order of the fit to the SS equation
if ~nargin
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    T_0 = result;
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_string      = @(x) ischar(x);

addRequired(  p, 'file_name',           val_string);
addRequired(  p, 'T_ref_string',        val_string);
addRequired(  p, 'T_string',            val_string);
addRequired(  p, 'SN',                  val_string);
addParamValue(p, 'profile_num',          default_profile_num,          val_numeric);
addParamValue(p, 'profile_min_P',        default_profile_min_P,        val_numeric);
addParamValue(p, 'profile_min_W',        default_profile_min_W,        val_numeric);
addParamValue(p, 'profile_min_duration', default_profile_min_duration, val_numeric);
addParamValue(p, 'order',        default_order,                        val_numeric);
addParamValue(p, 'vehicle_info', default_vehicle_info);

% Parse the arguments.
parse(p, file_name, T_ref_string, T_string, SN, varargin{:});

% Perform last stages of input validation.

profile_min_P        = p.Results.profile_min_P;
profile_min_W        = p.Results.profile_min_W;
profile_num          = p.Results.profile_num;
profile_min_duration = p.Results.profile_min_duration;
order                = p.Results.order;
vehicle_info         = p.Results.vehicle_info;

% end of input argument checking.
%%%%
% Get started by making sense of the inputs

% We need the name of the FP07 signals with and without pre-emphasis. The
% name without pre-emphasis is usually the section in the setup.cfg-file
% that contains the processing parameters. However, sometimes this signal
% does not exist and the information has to be gathered from the section
% for the signal with pre-emphasis. For eample, with a Sea-glider.

T_with_pre_emphasis_string = [T_string '_d' T_string]; % Should be T1_dT1 or T2_dT2
T_without_pre_emphasis_string = T_string; % Should be T1 or T2

warning off; % In case some variables do not exist
load(file_name, ...
    'setupfilestr', ...
    T_without_pre_emphasis_string, ...
    T_with_pre_emphasis_string, ...
    T_ref_string, ...
    't_fast', 't_slow', 'P_slow', 'W_slow', 'fs_fast', 'fs_slow')
warning on;

% We must identify the name of the section in the setup.cfg-file that
% contains the processing paramters.
section_name = T_without_pre_emphasis_string;

% In case that there is no thermistor signal without pre-emphasis
if exist(T_without_pre_emphasis_string,'var') 
    section_name = T_without_pre_emphasis_string;
    eval(['T_without_pre_emphasis = ' T_without_pre_emphasis_string ';']) 
else
    T_without_pre_emphasis_string = '[]';
    section_name = T_with_pre_emphasis_string;
    T_without_pre_emphasis = [];
end

% In case that there is no thermistor signal with pre-emphasis, we can only
% look at the signal without pre-emphasis and do not have to bother with
% deconvolution of the pre-emphasized signal.
if ~exist(T_with_pre_emphasis_string,'var') 
    T_with_pre_emphasis_string = '[]';
end

% make sure that there is at least one signal with FP07 data.
if ...
        isempty(T_without_pre_emphasis_string) && ...
        isempty(T_with_pre_emphasis_string)
    error(['Could not find any signals with base name = ' T_string])
end

% Finally, make sure that the temperature reference signal exists.
if ~exist(T_ref_string,'var')
    error(['Could not find any temperature reference signals with name = ' T_ref_string])
end


ratio = round(fs_fast / fs_slow); % sampling rate ratio
eval (['T_ref = ' T_ref_string ';']) % T_ref is the reference thermometer, usually SBT or JAC_T

eval(['T = ' T_without_pre_emphasis_string ';']) % T is the thermistor signal without pre-emphasis

% If we havce a signal with pre-emphasis, then we will use it to form the signal T.
if ~isempty(T_with_pre_emphasis_string)
    eval(['T = ' T_with_pre_emphasis_string ';']) % T is the thermistor signal with pre-emphasis
    T = deconvolve(...
        T_with_pre_emphasis_string, ...
        T_without_pre_emphasis, ...
        T, ...
        fs_fast, setupfilestr);
    T = reshape(T, ratio, []); % down size to match T_ref
    T = mean(T)';
end

title_string{1} = ['\rmFP07 \itin situ\rm calibration , SN-' SN];

% Figure out which section of the file to use for the in situ calibration.
if isempty(vehicle_info)
    profile_dir = 'down';
else
    profile_dir = vehicle_info.profile_dir;
end

    if exist('P_slow','var') && exist('W_slow','var')     
        if strcmpi(profile_dir, 'up') || strcmpi(profile_dir, 'down')
            % Get profile based on direction, duration, etc.
            profile = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, profile_dir, ...
                            profile_min_duration, fs_slow);                   
        elseif strcmpi(profile_dir,'glide')
            profile_down = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, 'down', ...
                            profile_min_duration, fs_slow);
            profile_up   = get_profile(P_slow, W_slow, profile_min_P, ...
                            profile_min_W, 'up',   ...
                            profile_min_duration, fs_slow);
            % Sort columns in ascending order
            profile = sort([profile_down profile_up],2);
        end
    end

    
% profile = get_profile(P_slow, W_slow, profile_min_P, profile_min_speed, direction, profile_min_duration, fs_slow);

profile_start = profile(1,profile_num);
profile_end   = profile(2,profile_num);
m = (profile_start:profile_end)';

%----------------------------------------------------
% --  Plot range for calibration -------------------- 
%----------------------------------------------------
figure(1), clf
plot(t_slow, T, t_slow(m), T(m),'r');grid on 
legend(...
     T_string, ...
    [T_string ' Profile'], 'location', 'northeast')
title (title_string)
xlabel('\it t \rm [s]')
ylabel('[counts]')

%-----------------------------------------------------
% -- Plot lag and corr between T  and T_ref ----------
%-----------------------------------------------------
figure(2), clf
junk_T      = detrend(T(m));
junk_T_ref  = detrend(T_ref(m));
range_T_ref = max(junk_T_ref) - min(junk_T_ref);
range_T     = max(junk_T) - min(junk_T);
junk_T      = junk_T * range_T_ref / range_T; % T should now span the same range as T_ref.

title_string{2} = [...
    'Detrended ' texstr(T_ref_string) ' & scaled ' texstr(T_string)];

plot(t_slow(m), [junk_T  junk_T_ref]);grid on
title (title_string)
xlabel('\it t \rm [s]')
ylabel('[ ^{\circ}C ]')
legend(texstr(T_string), texstr(T_ref_string))

%-----------------------------------------------------
% -- Plot xcorr coefficient for T and T_ref ----------
%-----------------------------------------------------
figure(3), clf
% Low-pass filtering the thermistor data to make it more compatible with the JAC-T
if strcmp(T_ref_string,'JAC_T')
    W_mean = abs(mean(W_slow(m)));
    fc = 0.73 * sqrt(W_mean / 0.62); % in Hz
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
else
    fc = fs_slow/3; % It is a Sea-Bird Thermometer
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
end

title_string{2} = texstr([...
    'X-correlation of ' T_string ...
    ' and ' T_ref_string]);
max_lag = round(10*fs_slow); % my estimate of the max lag required to find the actual lag.
[bb, aa] = butter(2,4/(fs_slow/2)); % 4 Hz smoother to suppress high-frequency noise

[correlation, lags] = xcorr(...
    filter(bb,aa,detrend(diff(T(m)))),...
    filter(bb,aa,detrend(diff(T_ref(m)))),max_lag,'coeff');
[max_corr, m_lag] = max(abs(correlation));
junk_m = m_lag; % needed for figure
m_lag = m_lag - max_lag - 1;

plot(lags/fs_slow, correlation, m_lag/fs_slow, correlation(junk_m), 'r*');grid on
xlabel('Lag [ s ]')
legend_string_1 = ['max X_{corr} = ' num2str(max_corr,2)];
legend_string_2 = ['@ \tau = ' num2str(m_lag/fs_slow,2) ' s'];
legend(legend_string_1, legend_string_2, 'location', 'northeast')
title(title_string)

%-----------------------------------------------------
% -- Do regression to get thermistor coefficients ----
%-----------------------------------------------------

% First align the T and T_ref signals using m_lag. 
% m_lag is expected to be negative.
T_ref2 = T_ref(m); % Only the profile and use a copy
T2     = T(m);

if m_lag >0, m_lag = 0; end

T_ref2 = T_ref2(1-m_lag:end);
T2     = T2(1:end+m_lag);
T_ref_regress = T_ref2 + 273.15; % in kelvin
T_ref_regress = 1 ./ T_ref_regress;

% Now gather information about the electronics for this thermistor.
my_object = setupstr( setupfilestr );
therm_type =    (char(setupstr( my_object, section_name, 'type')));
E_B = str2double(char(setupstr( my_object, section_name, 'E_B')));
a   = str2double(char(setupstr( my_object, section_name, 'a'  )));
b   = str2double(char(setupstr( my_object, section_name, 'b'  )));
G   = str2double(char(setupstr( my_object, section_name, 'G'  )));
adc_fs   = str2double(char(setupstr( my_object, section_name, 'adc_fs'  )));
adc_bits = str2double(char(setupstr( my_object, section_name, 'adc_bits'  )));
try zero = str2double(char(setupstr( my_object, section_name, 'adc_zero'  )));catch, zero = 0; end

if strcmp(therm_type, 'therm')
    factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);
    Z = factor*(T2 - a)/b;
elseif strcmp(therm_type, 't_ms')
    Z = T2 * (adc_fs/2^adc_bits) + zero; % Turn input into a voltage
    Z = ((Z - a)/b) *2 / (G*E_B);
end

%factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);

%Z = factor*(T2 - a)/b;
RT_R0 = (1 - Z) ./ (1 + Z); % This is the resistance ratio for this thermistor.
RT_R0 = log(RT_R0);

% Next confirm that the thermistor follows the Steinhart-Stein equation. The
%   plot should be a nearly straight line.

beta = zeros(1,order);
% Generate the coefficients for this thermistor.
p = polyfit(RT_R0, T_ref_regress, order);
pp = p; % save for later usage
p = 1 ./ p;
p = fliplr(p); % place in ascending order
T_0    = p(1);
for index = 2:order+1
    beta(index-1) = p(index);
end
Lag    = m_lag / fs_slow; % in seconds and should be negative.

% make a smooth line for comparison

R = linspace(min(RT_R0), max(RT_R0), 1000);
R = R';
T_inverse_predicted = polyval(pp,R);


%-----------------------------------------------------
% -- Plot regression ----------------------------------
%-----------------------------------------------------
figure(4), clf

h = plot(...
    T_ref_regress,       RT_R0, '.', ...
    T_inverse_predicted, R,     'r');
xlabel ('\itT \rm^{-1} [K^{-1}]')
ylabel ('log_{\ite\rm} (\itR_T\rm / \itR_{\rm0} \rm)')
set(h(1), 'markersize', 15)
set(h(2), 'linewidth',   2)

legend('Observed', 'Predicted','location', 'southeast')

title_string{2} = ...
    ['\itT_{\rm0}\rm = ' num2str(T_0) ', \beta_{\itn}\rm = ' num2str(beta) ];
title(title_string)
x_limits = get(gca,'xlim');
x_text = x_limits(1) + (x_limits(2) - x_limits(1))/25;
y_limits = get(gca,'ylim');
y_text = y_limits(2) - (y_limits(2) - y_limits(1))/10;
my_text = [...
    '$$\frac{1}{T} = \frac{1}{T_0} + ' ...
    '\frac{1}{\beta_1} \log_{\ e}\left(\frac{R_T}{R_0}\right) + ' ...
    '\frac{1}{\beta_2} \log^{\ 2}_{\ e} \left(\frac{R_T}{R_0}\right)$$'];

text(x_text, y_text, my_text, 'interpreter','latex', 'fontsize', 16)

%-----------------------------------------------------
% -- Use the computed co-efficients in a plot -------
%-----------------------------------------------------
figure(5), clf
fig_aspectratio(gcf,1.8);

if strcmp(therm_type, 'therm')
    factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);
    Z = factor*(T - a)/b;
elseif strcmp(therm_type, 't_ms')
    Z = T * (adc_fs/2^adc_bits) + zero; % Turn input into a voltage
    Z = ((Z - a)/b) *2 / (G*E_B);
end


%Z = factor*(T - a)/b;
RT_R0 = (1 - Z) ./ (1 + Z); % This is the resistance ratio for this thermistor.
RT_R0 = log(RT_R0);
pp = fliplr(1./p);

T_calibrated = polyval(pp, RT_R0);
T_calibrated = 1 ./ T_calibrated;
T_calibrated = T_calibrated - 273.15;

ax(1) = subplot(1,2,1);
plot(...
    T_ref(m-m_lag),  P_slow(m), ...
    T_calibrated(m), P_slow(m));grid on
set(gca,'ydir','rev')
xlabel('[ ^{\circ}C]')
ylabel('\itP   \rm[dBar]')
title(title_string{1})
legend(...
    texstr(T_ref_string), ...
    texstr(T_string), ...
    'location','southeast')

ax(2)=subplot(1,2,2);
c = get(0,'defaultaxescolororder');
plot(T_calibrated(m) - T_ref(m-m_lag), P_slow(m),'color',c(4,:));grid on
set(gca,'ydir','rev')
xlimits = get(gca,'xlim');
xlim(max(abs(xlimits))*[-1 1])
xlabel('[ ^{\circ}C]')
legend(texstr([T_string  ' - ' T_ref_string]),...
    'location','southeast')
% title(['\rm \itf_c\rm = ' num2str(fc,2) ' Hz'])
title('\rm Difference')
linkaxes(ax,'y')
 
 
 
 
