%% gradT_noise_odas
% Create a spectrum of the noise in a temperature gradient segment of data,
%  in physical units of (K/m)^2 / Hz.
%%
% <latex>\index{Functions!gradT\_noise\_odas}</latex>
%
%%% Syntax
%   gradT_noise = gradT_noise_odas (...
%           T_dT, T_dT_name, speed, F, setupfilestr, varargin, ...)
%
% * [T_dT]  A segemnt of raw pre-emphasized temperature data.
% * [T_dT_name] A string containing the name of the raw pre-emphasized data.
% * [speed] The speed of profiling during the time of the data segment, in m/s.
% * [F] A vector of frequency at which to evaluate the noise spectrum.
% * [setupfilestr] The string that contains the setupfile string for the
%      data file that contains the data segment. It is a standard product
%      of the conversion of a raw data file into a mat-file. Its name is
%      literally setupfilestr.
% * [noise_info] Optional structure containing the parameters of the noise
%      model for the raw data. If empty, default values are used.
% * [...] Optional configuration parameters to supplement, or override,
%           those values included within noise_info.
% * []
% * [gradT_noise] The noise spectrum for the identified data, in units of
%      (K/m)^2 / Hz.
%
%%% Description
% This function provides the frequency spectrum of the noise in a
% measurement of the fluctuations of the temperature gradient made with an
% FP07 thermistor mounted on an RSI instrument. The noise is estimated with
% the model described in Technical Note 040, Noise in Temperature Gradient
% Measurements.
%
% The default values of the parameters of the noise
% model are returned by calling this function without input arguments. They
% can be changed by presenting pairs of strings and values, individually.
% Or, they can be changed by placing the appropriate values into the fields
% of the structure $\texttt{noise\_info}$.
%
% $\texttt{T\_dT}$ is a segment of data for which you want the noise
% spectrum. The segment is used to estimate the raw to physical conversion
% factor, only. The string that identifies the channel of the
% pre-emphasized data is usually called $\texttt{T1\_dT1}$ or
% $\texttt{T2\_dT2}$. You can also use the abreviated names of
% $\texttt{T1}$ and $\texttt{T1}$. The name is needed to find the
% electronic and sensor parameters for the data segment.
%
% $\texttt{speed}$ can be a vector of length equal to the length of
% $\texttt{T\_dT}$, or a scalar. If it is a vector, the average value is
% used to estimate a single value for speed. $\texttt{F}$ is a vector of
% frequency at which to estimate the spectrum of the temperature-gradient
% noise.
%
% $\texttt{setupfilestr}$ is the string that contains the setup information
% for the data file that contains the data segment. It is a standard
% product of the convserion of a raw data file (a $\texttt{*.p}$ file) into
% a $\texttt{*.mat}$ file in which the data are in physical units.
%
% If $\texttt{gradT\_noise\_odas}$ is called without input arguments, it
% returns the default parameters used in the noise model. For example,
%
%    >> noise_info = gradT_noise_odas
%
%    noise_info =
%               a_0: 1.4e-17
%               f_a: 3
%               b_0: 5.6e-17
%               f_b: 47
%                 K: 1.38e-23
%                 R: 3000
%              f_AA: 110
%       make-figure: false
%
% if make_figure is true, then the noise spectrum is plotted.
%

% * 2017-05-05 RGL Original
% * 2017-06-08 RGL Brought to functioning form.
% * 2018-03-08 JMM Fixed a bug so name_without_preemphasis is correct



function gradT_noise = ...
    gradT_noise_odas (T_dT, T_dT_name, speed, F, setupfilestr, varargin)


% Default values used to model the noise of a pre-emphasized temperature
% channel
default_a_0         = 1.4e-17; % stage 1 , HF noise
default_f_a         = 3; % stage 1 knee frequency
default_b_0         = 5.6e-17; % stage 2 HF noise
default_f_b         = 47; % stage 2 knee frequency
default_K           = 1.38e-23; % Boltzman's constant
default_R_0         = 3000; % Bridge resistors
default_f_AA        = 110; % Anti-aliasing filter parameter.
default_make_figure = false; % do not plot the spectrum

if ~nargin
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        gradT_noise.(char(param{1})) = eval(d.name);
    end
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x)  isnumeric(x) && isscalar(x) &&    (x) > 0;
val_speed       = @(x)  isnumeric(x)                && all(x) > 0;
val_vector      = @(x) (isnumeric(x) && isvector(x));
val_string      = @(x)  ischar   (x);
val_logic       = @(x)  islogical(x);

addRequired(  p, 'T_dT',                               val_vector);
addRequired(  p, 'T_dT_name',                          val_string);
addRequired(  p, 'speed',                              val_speed);
addRequired(  p, 'F',                                  val_vector);
addRequired(  p, 'setupfilestr',                       val_string);
addParamValue(p, 'a_0',           default_a_0,         val_numeric);
addParamValue(p, 'f_a',           default_f_a,         val_numeric);
addParamValue(p, 'b_0',           default_b_0,         val_numeric);
addParamValue(p, 'f_b',           default_f_b,         val_numeric);
addParamValue(p, 'K',             default_K,           val_numeric);
addParamValue(p, 'R_0',           default_R_0,         val_numeric);
addParamValue(p, 'f_AA',          default_f_AA,        val_numeric);
addParamValue(p, 'make_figure',   default_make_figure, val_logic);

% Parse the arguments.
parse(p, T_dT, T_dT_name, speed, F, setupfilestr, varargin{:});


T_dT          = p.Results.T_dT;
T_dT_name     = p.Results.T_dT_name;
speed         = p.Results.speed;
F             = p.Results.F;
setupfilestr  = p.Results.setupfilestr;
a_0           = p.Results.a_0;
f_a           = p.Results.f_a;
b_0           = p.Results.b_0;
f_b           = p.Results.f_b;
K             = p.Results.K;
R_0           = p.Results.R_0;
f_AA          = p.Results.f_AA;
make_figure   = p.Results.make_figure;

speed = mean(speed); % Only need the mean speed

%---------------------------------------------------------
% Now decipher the name of the thermistor signal
if length(T_dT_name) <= 2 % we have a shortenned name
    name_with_pre_emphasis = [ T_dT_name '_d' T_dT_name];
    name_without_pre_emphasis = T_dT_name;
end

if length(T_dT_name) > 2 % we may have a name of a pre-emphasized channel
    index_to_underscore = strfind(T_dT_name, '_');
    if ~isempty(index_to_underscore)
        index_to_underscore = index_to_underscore(1)-1; % in case of more than one.
        name_with_pre_emphasis = T_dT_name;
        name_without_pre_emphasis = T_dT_name(1:index_to_underscore);
    else
        error(['Cannot decipher the name ' T_dT_name])
    end
end

obj = setupstr(setupfilestr); % make object for increased speed

if isempty(setupstr(obj, '', 'name', name_without_pre_emphasis))
    name_without_pre_emphasis = [];
end

if isempty(setupstr(obj, '', 'name', name_with_pre_emphasis))
    error(['A channel with name ' ...
        name_with_pre_emphasis ...
        ' is not in the configuration string'])
end


%---------------------------------------------------------
% Now lets get the basic parameters that we need
fs = setupstr( obj, 'root','rate');
if ~isempty(fs)
    fs = str2double(fs);
else
    error('Cannot extract sampling rate, fs, out of the setupfile string')
end

diff_gain = setupstr( obj, name_with_pre_emphasis, 'diff_gain' );
if isempty(diff_gain)
    error(['There is no diff_gain parameter for the channel ' ...
        name_with_pre_emphasis])
else
    diff_gain = str2double(diff_gain);
end

% Make a cell array of the parameters that we need and the optional ones.
params_we_need = {...
    'a', 'b', 'adc_fs', 'adc_bits', 'g', 'e_b', 't_0', ...
    'beta', 'beta_1', 'beta_2', 'zero', 'type'};

beta = {[]}; % This is needed because beta is a Matlab built-in function.
type = {[]}; % This is needed because type is a Matlab built-in function.

% First check the section of the setup-file string for the channel with
% pre-emphasis.
name = name_with_pre_emphasis; % make it shorter
for index = 1 : length(params_we_need)
     junk = setupstr( obj, name, params_we_need(index));
     eval([params_we_need{index} ' = junk;'])
end

%---------------------------------------------------------
% Now go to the section without pre-emphsis to get the parameters that may
% still be missing.
if ~isempty (name_without_pre_emphasis)
    name = name_without_pre_emphasis; % make it shorter
    for index = 1 : length(params_we_need)
        if isempty (eval(params_we_need{index}))
            junk = setupstr( obj, name, params_we_need(index));
            eval([params_we_need{index} ' = junk;'])
        end
    end
end

if ~isempty(beta) && isempty(beta_1)
    beta_1 = beta;
end

if isempty (zero)
    zero = 0;
else
    zero = str2double(zero);
end
if isempty (beta_2)
    beta_2 = inf;
else
    beta_2 = str2double(beta_2);
end

a        = str2double(a);
b        = str2double(b);
adc_fs   = str2double(adc_fs);
adc_bits = str2double(adc_bits);
g        = str2double(g);
e_b      = str2double(e_b);
T_0      = str2double(t_0);
beta_1   = str2double(beta_1);


%---------------------------------------------------------
% deconvolve T_dT to remove pre-emphasis

T = deconvolve (name_with_pre_emphasis, [], T_dT, fs, obj);

% Now compute the absolute temperature
if strcmpi(type, 'therm')
    Z = ((T -a) / b) * (adc_fs / pow2(adc_bits)) * 2 / (g * e_b);
end
if strcmpi(type, 't_ms')
    Z = T * (adc_fs / pow2(adc_bits)) + zero;
    Z = ((Z - a) / b) *2 / (g * e_b);
end

R = (1 - Z) ./ (1 + Z); % resistance ratio, R_T/R_0
R(R<0.1) = 1; % in case of broken thermistor
log_R = log(R);

T = 1/T_0 + log_R / beta_1;
if ~isfinite(beta_2)
    T = T + log_R.^2 / beta_2;
end
T = 1./ T; % now in Kelvin
T = mean(T);
T_celsius = T - 273.15;

r = mean(R);

%---------------------------------------------------------
% Calculate the scale factor for the conversion of the temperature
% gradient.

    eta = (b / 2) * pow2(adc_bits) * g * e_b / adc_fs;
    scale_factor = 1;
    if isfinite(beta_2)
        scale_factor = 1 + 2 * (beta_1 / beta_2) * log_R;
    end
    scale_factor = scale_factor .* T.^2 .* (1 + R).^2 ./ (2 * eta * beta_1 * R);

    scale_factor = scale_factor ./ speed; % Convert the time-derivative to a gradient.
    scale_factor = mean(scale_factor);

%---------------------------------------------------------
% Now calculate the thermistor noise model

G_D     = diff_gain;
G_1     = g;
delta_s = adc_fs / pow2(adc_bits);
R       = mean(R * R_0);
phi_R   = 4 * K * T .* R;

Noise_1     = G_1^2 * (2*a_0 * (1 + f_a ./ F) + phi_R);

G_2         = 1 + (2 * pi * G_D * F).^2;

Noise_2     = G_2 .* (Noise_1 + 2*b_0 * (1 + f_b ./ F));

G_AA        = (1 + (F/f_AA).^8 ).^2; G_AA = 1 ./ G_AA;

Noise_3     = Noise_2 .* G_AA;

Noise_4     = Noise_3 + 2.7 * (2 * delta_s^2) ./ (12 * fs);

Noise_N     = Noise_4 / delta_s^2;

G_HP        = (2*pi*G_D*F).^2 ./ (1 + (2*pi*G_D*F).^2);

gradN_noise = Noise_N .* G_HP;

gradT_noise = gradN_noise * scale_factor^2;

%---------------------------------------------------------
% plot results, if requested
if make_figure
figure(1000)
clf
x_lim = [0 0];
if F(1) == 0
    x_lim(1) = F(2);
else
    x_lim(1) = F(1);
end
x_lim(2) = F(end);


h = loglog(F, gradT_noise, 'linewidth', 2);
set(h(1),'color','k')
xlabel('\itf\rm  [Hz]')
ylabel('[ K^2 m^{-2} Hz^{-1}]')
legend ('Thermistor Noise \nabla\itT_{\rm }', 'location', 'northwest')
set(gca, 'xlim', x_lim)

title_string = { ...
    ['\rmRSI Thermistor Noise, ' name_without_pre_emphasis], ...
    ['\itT\rm = ' num2str(T_celsius,4) '\circC, \itR_T\rm /\itR_{\rm0}\rm = ' num2str(r,3) ...
    ', ' num2str(speed,3) ' m s^{-1}']};

title(title_string)
end
