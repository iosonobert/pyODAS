%% odas_version_info - EXCLUDED
% Return the version of the current ODAS library.
%%
% <latex>\index{Functions!odas\_version\_info}</latex>
%
%%% Syntax
%   version = odas\_version\_info( )
%
% * []
% * [version] Numeric value of the version.
%
%%% Description
% Function called by other functions to determin the current version of
% the ODAS Libarary.  Used when checking for outdated data files.

% Version History
%
% * 2015-06-03 (WID) initial version 4.0
% * 2016-08-02 (WID) updated to version 4.1.  Added support for some JAC
%                    sensors.  Minor bug fixes.  Notable inclusion of
%                    show_ch fuction.
% * 2017-01-27 (WID) updated to version 4.2.  
%                    addition of AEM1-G analog and digital types
% * 2018-01-05 (JMM) updated to version 4.3.
%                    Bug fixes, including despiking of piezo-accelerometers
%                    This version will be the last major update to the ODAS
%                    library. Minor improvements and bug fixes will be
%                    made, but major changes will be released as part of
%                    the Zissou software package.
% * 2018-04-14 (JMM) updated to version 4.3.02
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at Lake Tahoe).
% * 2018-04-25 (JMM) updated to version 4.3.03
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at UNAL). Most significant
%                    change is that a plot of pressure will be generated
%                    even if no profiles are detected. 

function version = odas_version_info( )

version = 4.3;

end

