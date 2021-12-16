%% show_ch
% Extract approximate channel and time vectors from a RSI raw binary data
% file.
%%
% <latex>\index{Functions!show\_ch}</latex>
%
%%% Syntax
%
%   [ch, t] = show_ch( filename, channel, range, convert )
%
% * [filename] - String containing the name of a RSI raw binary data file.
%                When undeclared or left empty this function will prompt
%                for a file name. The default value provided at this prompt
%                is the most recently modified data file in the current
%                directory.
% * [channel] - Channel name as specified within the configuration file
%               used when collecting this data file.  Default value is "P".
% * [range  ] - value pair, [start end] for the first and last record to be
%               returned or plotted.
% * [convert] - Logical value that determines if data are converted into
%               physical units.  Default is "true" or do convert.
% * []
% * [ch] - Data or vector of the specified channel. When possible, this 
%          data will be converted into physical units via the deconvolve 
%          and convert_odas functions.
% * [t] -  Time vector associated with the specified channel.
%
%%% Description
% Extract a specific channel from an RSI raw binary data file.  Both the
% channel data and associated time vector are returned.  By default, a plot
% of the returned data / time vector is generated allowing for a quick
% preview of channel data within a data file.
%
% Resulting data vectors are automatically deconvolved and converted into
% physical units. This function notifies the user when these conversions
% are automatically applied.
%
% It should be noted that this function bypasses any previously generated
% .MAT files and accesses the RSI raw binary file directly. No
% modifications are made when this function is called making it a safe
% function to use when quickly examining a data file.
%
% This function only works on RSI raw binary data files version 6 or
% greater.
%
%%% Examples
%
%     >> show_ch;
%
% Plot the pressure channel of a specified data file. Because a data file
% was not specified directly, this function prompts for a file name.
%
%     >> [T1,t_T1] = show_ch( 'DAT_001.P', 'T1' );
%
% Extract and convert the channel named "T1" from the data file
% "DAT_001.P". Note the resuting figure has been supressed with values of 
% T1 and t_T1 being returned in place of the figure.
%
%     >> show_ch( 'DAT_001.P', 'T*' );
%
% Plot all channels that start with the letter "T" from the data file 
% "DAT_001.P".
%
%     >> show_ch( 'DAT_001', {'T1', 'sh1'}, false );
%
% Extract and plot channels "T1" and "sh1" from the data file "DAT_001". 
% Note that the resulting figure has data vectors plotted as raw counts.
%

% *Version History:*
%
% * 2016-06-27 (WID) Initial version.
% * 2016-08-02 (WID) Empty channel name default to "P"
% * 2016-10-06 (WID) Instigated the option of no plotting.
% * 2016-12-09 RGL, added input argument for no conversion into physical
%              units.
% * 2017-01-27 (WID) Removed the plot option. Now either plot or return
%              values but not both.  When reading very large data files,
%              extract the data portion before processing.
% * 2017-05-09 (WID) Default range starts at 0 - previously set to 1.
% * 2018-01-18 (WID) Better support instruments having a channel with
%              pre-emphasis but no channel without pre-emphasis.

function [a, b] = show_ch( fname, channel_names, range, convert)

if nargin < 3, range = []; convert = true; end
if nargin < 4, convert = true; end
do_plot = nargout == 0;

% Set default range to be the entire file.
if isempty(range), range = [0,Inf]; end

tempfilename = '';
padded_time = 0;

if nargin < 1 || isempty(fname)
    [vars, d] = read_odas();
else
    
    workingfile = fname;
    
    % Check to see of the range is such that we benefit from creating a
    % temporary file.  If so, create the file and adjust the range to work
    % with the new file.
    if range(1) > 60 || range(end) < Inf
        extract_start = max(range(1) - 60, 1);
        extract_stop  = range(end);
        padding_start = range(1) - extract_start;
        range = [padding_start, range(2) - range(1) + padding_start];
        tempfilename = [tempname, '.p'];
        extract_odas( fname, extract_start, extract_stop, tempfilename );
        workingfile = tempfilename;
        padded_time = extract_start;
    end
    
    [vars, d] = read_odas( workingfile );
end

% Remove the temporary file - but only if we created one.
if ~isempty(tempfilename)
    delete(tempfilename);
end

% Ensure the input channel name(s) are formatted as a cellstr array.
if nargin < 2 || isempty(channel_names)
    channel_names = {'.*'};
elseif ischar(channel_names)
    channel_names = {channel_names};
elseif ~iscellstr(channel_names)
    error('Channel name must be either a string or a cell array.');
end

% Expand wildchar characters to select all maching channel names.
channels = {};
for cname = channel_names
    sections = setupstr(d.cfgobj, cname);
    if isempty(sections)
        channels(end+1) = cname;
    else
        for section = sections
            fixed_name = setupstr(d.cfgobj, section, 'name');
            cmp_result = regexp(char(fixed_name), char(cname));
            if ~isempty(fixed_name) && (isempty(cmp_result) || (cmp_result == 1))
                channels(end+1) = fixed_name;
            end
        end
    end
end


if isempty(channels)
    error('Unable to find channel name: %s within the data file.', channel_names{1});
end


% The first plot sets the figure / axes settings.  Additional plots are
% much simpler.
legend_str = {};

previous_channels = {};

% If not plotting, extract the first matching channel name then exit.
[ch,t,units,chname] = extract_ch( channels{1}, previous_channels, convert);
n = find( t>=range(1) & t<=range(end) );

% When exiting, return the data if requested.  If two output vectors are
% requested, arrange the time/data vectors so they are compatible with the
% plot command.
if ~do_plot
    if nargout >= 1, a = ch(n); end
    if nargout >= 2, b = t(n); end
    return;
end

if ~isempty(chname), previous_channels(end+1) = {chname}; end

defer_pressure = [];
figure();
h = axes();
grid('on');
xlabel('Time [s]', 'Interpreter', 'tex');

% Only add a vertical axis label if there is only 1 channel to plot.
if length(channels) == 1
    
    plot(h, t(n) + padded_time, ch(n));
    legend_str(end + 1) = {texstr(chname)};
    
    if strcmpi(chname, 'P')
        % Reverse the y-axis if this is a pressure plot.
        set(gca, 'ydir', 'rev');
        % Type poly does not have units by default.
        if isempty(units) || length(units) == 1, units = {'[dBar]'}; end
    end
    
    % If no units were found, set to an empty string.
    if isempty(units), units = {''}; end
    ylabel([texstr(chname) ' ' char(units)], 'Interpreter', 'tex');
    
else
    
    hold('on');
    label_units = 'first pass';
    for c = 1:length(channels) %channels(1:end)
        if c == 1
            time = t;
            dat = ch;
        else
            [dat,time,units,chname] = extract_ch( channels{c}, previous_channels, convert);
            if isempty(chname)
                continue;
            else
                previous_channels(end+1) = {chname};
            end
        end
        
        n = find( time>=range(1) & time<=range(end) );
        
        % We print units for the vertical axis only if all channels have
        % the same units.
        if strcmpi('first pass', label_units) || strcmpi(units, label_units)
            label_units = units;
        else
            label_units = '';
        end
        
        if strcmpi(chname, 'P')
            defer_pressure.dat = dat(n);
            defer_pressure.time = time(n) + padded_time;
            defer_pressure.units = units;
            defer_pressure.chname = chname;
        elseif ~any(strcmp(legend_str, chname))
            legend_str(end + 1) = {texstr(chname)};
            plot(h, time(n) + padded_time, dat(n));
        end
    end
    
    ylabel(label_units, 'Interpreter', 'tex');
    
    % Plot pressure last, so it overlaps the other channels.
    if ~isempty(defer_pressure)
        legend_str(end + 1) = {texstr(defer_pressure.chname)};
        plot(h, defer_pressure.time, defer_pressure.dat, 'LineWidth', 2, ...
            'Color', 'black');
    end
    
    hold('off');
    
end

% Add a brief title.
[P,N,E] = fileparts(d.fullPath);
title(gca, texstr(['show_ch plot: ' N E]));

legend(legend_str);





    function [ch, t, units, n] = extract_ch( channel_name, previous_channels, convert)
        
        if ~isfield(d, channel_name)
            error('Channel "%s" not found within data file "%s".', ...
                channel_name, fname);
        end
        
        
        %%%% Perform a deconvole if required.  It is required when the name matches
        %%%% the "X_dX' naming convention and a 'diff_gain' parameter is found.
        n = channel_name;
        L = length(n);
        if L >= 4 && mod(L,2) == 0 && ...
                strcmp(n(L/2:L/2+1), '_d') && ...
                strcmp(n(1:L/2-1), n(L/2+2:end))
            ndn = n;
            n = n(1:L/2-1);
        else
            ndn = [n '_d' n];
        end
        
        data = [];
        if isfield(d, n), data = d.(n); end
        
        
        % This channel has already been extracted. Do not do it again.
        if any(strcmp(n,previous_channels))
            ch = [];
            t = [];
            units = '';
            n = '';
            return;
        end
        
        
        %%%% When determining the length of the data vector, use the longest value
        %%%% (with/without preemphasis).
        chlength = 0;
        if isfield(d, ndn)
            chlength = length(d.(ndn));
        elseif isfield(d, n)
            chlength = max(length(d.(n)), chlength);
        else
            error('Unable to find the specified channel: %s', channel_name);
        end
        
        % Generate / extract the time vector.
        if chlength == length(d.t_fast)
            t = d.t_fast;
        elseif chlength == length(d.t_slow)
            t = d.t_slow;
        else
            % Neither a fast nor slow channel.  Calculate a new time vector.  Note
            % that it will be slightly off - a fraction of a second for the entire
            % profile.  No one will ever notice...
            t = (1:length(d.(channel_name)))-1;
            t = t * d.t_slow(end) / length(d.(channel_name));
            warning('Time vector a close approximation, not an exact value.');
        end
        
        
        %%%% Go back and finish off the deconvolve.  Previously, we never had the
        %%%% "fs" argument so processing was delayed.
        try
            fs = length(t) / length(d.t_fast) * d.fs_fast;
            data = deconvolve(ndn, data, d.(ndn), fs, d.cfgobj);
            disp(['Performing a deconvolve on channel ' ndn]);
        catch
            data = d.(channel_name);
            n = channel_name;
        end
        
        
        %%%% Convert the data vector into physical units.
        ch = data;
        units = '[counts]';
        if convert
            try
                [ch, units] = convert_odas( data, {ndn,n}, d.cfgobj );
                disp(['Converting channel ' n ' into physical units.']);
            catch
                warning('Unable to convert channel "%s".  Returning units of [counts].', n);
            end
        end
        
        %%%% Are units empty?  If so, attempt to get them from the configuration
        %%%% file.
        if strcmp(units, ' ')
            units = char(setupstr(d.cfgobj, n, 'units'));
        end
        
    end


end