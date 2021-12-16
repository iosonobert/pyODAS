%% plot_HMP
% Real-time or replay plotting of data collected with a horizontal microstruture
% profiler (HMP).
%%
% <latex>\index{Functions!plot\_HMP}</latex>
%
%%% Syntax
%   plot_HMP( fileName, setup_fn )
%
% * [fileName] name of the data file
% * [setup_fn] the name of the channel/data setup file
%
%%% WARNING
% This function is very old and does not make use of modern data files.  It
% is included for users who already use the function and is not recommended
% for new users.
%
%%% Description
% Function for the real-time (or post-time) plotting of data collected with
% a horizontal profiler. The plotting style is in the form of horizontal
% traces in a stack of subplots as shown in the next figure. Up to 16
% channels can be plotted. Each screen consists of 50 records of data
% (which are usually 50 seconds long).  All graphical data are displayed in
% units of counts. The channels that are plotted and the names placed on
% the y-axis labels are determined by entries in the configuration file. In
% addition to the time-series, the function will also print on to the plot
% the average temperature from a Sea-Bird thermometer, the pressure, tow
% speed, and heading, if these are available. The averages are for the
% first record on the figure. The name of the data file, the date and time
% when the data were collected and when they were plotted are also
% annotated into the figure. Plotting continues for as long as new data are
% available in the source file. The function is somewhat interactive and
% the user can choose to have the local printer produce a hardcopy of each
% 50-record figure. The user can also select to look at a particular
% 50-record segment. 
%
%%% Examples
%
% The next Figure was produced with the following parameters located in the
% setup file:
%
%     plotting: 1,2,3,7,49,50,51,10,11,53,54,14,56,55,16,17,18,19
%     plotnames: Ax,Ay,Az,T2\_dT2,F_R,F_N,F_C,P,P\_dP,U,V,T7\_dT7,Alt_error,
%                Alt,SBT2E,SBT2O,SBC2E,SBC2O
%     plotaverages: 16,17,10,53,54,32,33
%     averagenames: SBT2E,SBT2O,Pres,U,V,Mx,My
%
% @image @images/plot_HMP @Example output produced by $\texttt{plot\_HMP}$.
% @Example output produced by $\texttt{plot\_HMP}$. Data courtesy of Eric
% Kunze. 
%
%   Note: This function has not been fully converted to support ODAS v6 and 
%   higher. In the meantime, please note the following guidelines and 
%   limitations (only for ODAS v6 and higher):
%
%   * There is no need to supply a setup file name. All plotting related
%     parameters are hardcoded at the top of this function or obtained
%     from the setup file string.
%   * Channel averaging is only supported for Pressure at this time.
%   * Channel averaging for Pressure is configured by hardcoding three 
%     variables (found at the top of this m-file) as follows:
%         ch_mean_seq = [10]; 
%         ch_mean_names = {'Pres'};
%         ch_mean_sections = {'pres'};
%   * Plotting of channels is accomplished by hardcoding the variable 
%     ch_plot_seq at the top of this function, for example:
%         ch_plot_seq = [1 2 3 4 5 6 7 8 9 10 11]';
%

% *Version History:*
%
% * 1994-11-10 (LZ) initial
% * 1995-02-01 (RAMO) use static read_block10.m file
% * 1998-03-01 (FW) major revisions
% * 1999-06-01 (RGL) Minor changes
% * 2004-06-01 (IG) Existing parts of the routine separated in M-files
% input_blocks.m, plot_ch_data.m,calibrations_plotting.m, load_setup_data.m.
% Now read calibration parameters directly from ODAS setup file.  Display is
% updated every ten seconds.  Added sections of code throughout sub-routines to
% optionally deal with the new scount board.  Other comments regarding changes
% are inserted in the code below and in the sub-routines.
% * 2007-09-01 (RGL) Added Alec EM Current meter, and PNI Magnetomer
% * 2009-03-06 (RGL) use fopen_odas instead of fopen
% * 2011-05-18 (AWS) Moved the following m-files back into plot_HMP:
% assign_plot_ch.m, calibrations_plotting.m, input_blocks.m, load_setup_data.m,
% plot_ch_data.m, read_block10.m. These functions are only used by plot_HMP.
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-04-11 (WID) replaced inifile_with_instring calls with setupstr
% * 2012-09-09 (WID) updated documentation to allow for publishing
% * 2012-11-05 (WID) updated documentation
% * 2015-07-03 (WID) quick look of code showed some errors - corrected
%                    them.


function plot_HMP(fileName, setup_fn)
%**************************************************************************
%                                                                         *
%  Notes:                                                                 *
%     This function was written using Matlab and designed for real time   *
%     data acquisition. For the usage, please see the documentation for   *
%     details.                                                            *
%                                                                         *
%**************************************************************************

% parsing of command line parameters
if nargin < 2; 
    setup_fn   = 'setup.txt';
    disp('Using local setup.txt file')
end
if nargin < 1; fileName  = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization of constants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Defaultaxesfontsize',10);
set(0,'Defaulttextfontsize',10);

size_of_integer =  2;% size of integers from data files in bytes
header_size_i   = 18;% index of header_size in header record
block_size_i    = 19;% index of block_size in header record

header_in_byte = 128;% the size of the header of each block (in byte)
plot_blocks    =  50;% number of blocks to plot on screen each time 
pts_per_block  =  32;% number of data points of each block after decimation
pts_per_ch     = plot_blocks * pts_per_block;% total number of data points of each channel to plot
real_time      =   0;% set flag for real time, default = not real time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Input from the User
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real_time = 1;
answer = input('Is this real time plotting (n/y): (default = y) ', 's');
if strcmp(answer, 'n')
  real_time = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % NOTE:
  %     (nb = 0) means that plotting data will not
  %     stop until the end of file is reached
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nb = 0;
end

% some parts of start_odas added from here by RAMO
header_size = 64;% the size of the header of each block (in integers)

if isempty(fileName); fileName = getDataFileName; end; % get data file name from the user
if (fileName == 'q'); return;end

% load channel data file by calling 'load_ch_data' function
% IG, June 2004: now read the main ODAS setup file, rather than two
% separate parameter files specific to this routine
[ch_plot_seq,ch_names,ch_mean_seq,ch_mean_names,ch_pars,two_ch_var]=load_setup_data(setup_fn);
max_ch = length(ch_plot_seq);


% extract some record parameters    
fid = fileOpen(fileName);
HD = fread(fid, header_size, 'ushort');	% read header
frewind(fid);			% back to the beginning of the file
header_in_byte = HD(header_size_i);
block_size_byte = HD(block_size_i)-header_in_byte;

header_version = bitshift(HD(11), -8) + bitand(HD(11), 255) /1000;
if header_version >= 6
    first_record_size_in_bytes = header_in_byte + HD(12);
else
    first_record_size_in_bytes = block_size_byte + header_in_byte;
end


% compute sizes of setup matrix, total number of channels
[rows, cols, no_slow_cols, slow_ch, fast_ch] = load_ch_setup(fid);
no_of_ch = length(slow_ch) + length(fast_ch);	

max_rows = block_size_byte / (size_of_integer * cols);
records_per_block = max_rows / rows;

seq      = assign_plot_seq(ch_plot_seq, slow_ch, fast_ch); % assign the plotting channels
mean_seq = assign_plot_seq(ch_mean_seq, slow_ch, fast_ch); % assign the plotting channels

% RAMO added this part to replace writing out read_block10 each time

% paramv indices
header_size_i=1;        % index of the header size parameter 
rows_i =2;              %          number of rows of channels
cols_i =3;              %          number of columns of channels
num_slow_col_i=4;       %          number of slow channel columns
recs_per_block_i=5;     %          records per block
points_wanted_i=6;      %          points wanted (ie decimation factor)

paramv(header_size_i) = header_size;
paramv(rows_i) = rows;
paramv(cols_i) = cols;
paramv(num_slow_col_i)= no_slow_cols;
paramv(recs_per_block_i)= records_per_block;
paramv(points_wanted_i)= pts_per_block;

% end of RAMO change

% end of  parts from start_odas moved here by RAMO

[start_block,end_block,nb]=input_blocks(real_time,0,plot_blocks,0);

% set flag for printing this figure, default = no writing 
hardcopies = 0;
answer = input('Do you want hardcopies (n/y): (default = n) ', 's');
if strcmp(answer, 'y')
   hardcopies = 1; 
end

% set file pointer to the 'start_block'
if start_block >= 0
   block_size = ...
      rows * cols * records_per_block * size_of_integer + header_in_byte; 
   status = fseek(fid, first_record_size_in_bytes + (start_block)*block_size, 'bof');
   if status ~= 0
      error ('*** The value of starting block might be too large ***'); 
   end
end

H = zeros(1, header_size);% preallocate the buffer 
S = zeros(1,pts_per_ch*no_of_ch);
S = reshape(S, pts_per_ch, no_of_ch);

nn = 0;% current plotting block data number 
nrem = 0;% which block 10
done = 0;% flag for all done
read_success = 1;% flag for reading status
continued_flag=0; % whether or not are continuing previous plotting (determines how update plots)
skip_dist = 0;  % If have continued, how many blocks have we skipped

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Program Main Loop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
while ~done 

   status = 0; % set reading flag
   nrem = rem(nn,5);  % which block 10
   if nn == 0 && ~continued_flag; S = S * NaN;	end; % initialize S only at the first block 10

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %    Main Reading Loop
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % keep reading data if there is more data (status==0)
   % but stop reading  when 5 ten block channels are read 
   % i.e. when (nrem > 4) occurs

   while status == 0 && nrem <= 4
       if nn>=1 || continued_flag, %nn>4 || continued_flag,         
% Sliding window of points if trying to look at more than 50 points at a time
           if nrem==0 && continued_flag && skip_dist>0
               S(1:pts_per_ch-skip_dist,:) = S(skip_dist+1:end,:);
               S(pts_per_ch-skip_dist+1:end,:) = NaN;
           end
           S(1:pts_per_ch*4/5,:) = S(pts_per_ch*1/5+1:end,:);
       end
       [status, S, H] = read_block10(fid, S, H, 4, paramv);
  
      if status >= 0; % reading successful
         nn = nn + 1; % get ready to read the next block
         if (nb ~= 0 && ~done);% not real time and not done
            done = (nn >= nb);
         end

         if nrem == 4 || done; break; end;% got 5 blocks 10 or done, then exit reading loop
         if status == 0; nrem = nrem + 1; end; % (status == 0) means there is data in the data file


      % reading failed, means the end of file is reached
      % find out if there are more files to read
      elseif real_time 
         read_success = 0;		% reading failed
         nrem = nrem - 1;
         if nrem < 0
            fprintf('**** Warning: reached the end of data file %s. ****\n', fileName);
         end
      else
         nrem = nrem - 1 ;
         done = 1;
         read_success = 0;
            fprintf('**** Warning: reached the end of data file %s. ****\n', fileName);
      end % if status >=0
      
      if nrem >=0&&read_success~=0
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% plotting channel data on screen
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot_ch_data(fileName,S,H,nn,4,seq,ch_names,mean_seq,ch_mean_names,ch_pars,pts_per_block,pts_per_ch,...
            round(start_block-plot_blocks*(1-nrem/5)),plot_blocks,0,two_ch_var)
    end
end % while status == 0 && nrem <= 4
   
if nrem >= 0  && read_success~=0
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% plotting channel data on screen
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	day_flag = (nrem == 4 || done || ~read_success);
    plot_ch_data(fileName,S,H,nn,4,seq,ch_names,mean_seq,ch_mean_names,ch_pars,pts_per_block,pts_per_ch,...
        round(start_block-plot_blocks*(1-(nrem+1)/5)),plot_blocks,day_flag,two_ch_var);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	%   Do we send it to the plotter/printer?
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% only at 5th ten block segment or all done
      if nrem == 4 || done || ~read_success
         if hardcopies; %             print
              print(gcf,'-dps2', [fileName '_' num2str(start_block) '.ps'])         
         end 
     end
 end % nrem >= 0

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %     Move on to the next set of blocks
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
continued_flag=0;       % Whether or not to continue the plot on the next run
old_endblock=end_block;
   if ~real_time && done
       [start_block,end_block,nb,quit_flag]=input_blocks(real_time,done,plot_blocks,end_block);
       if ~quit_flag
         if start_block >= 0
            % skip over to the starting block
            block_size = rows * cols * records_per_block * 2 + header_in_byte; 
            status = fseek(fid, first_record_size_in_bytes + (start_block)*block_size, 'bof');
            if status ~= 0
               error ('***** The value of starting block might be too large *****'); 
            end
         end
         nn = 0;  % initialization variables
         nrem = 0;
         done = 0;
     else fclose(fid);
     end % strcmp(answer, 'y')
   elseif real_time && ~read_success
      fclose(fid);     % Close the existing file (added by IG, June 2004)
      fileName = getDataFileName;  % get a new file name and open new file
      if (fileName == 'q'); return;end
      fid = fileOpen(fileName);
      
      % open succeds, set flag
      [start_block,end_block,nb]=input_blocks(real_time,done,plot_blocks,0);
               
      % set file pointer to the 'start_block' 
      if start_block >= 0
         block_size = rows * cols * records_per_block * 2 + header_in_byte; 
         status = fseek(fid, first_record_size_in_bytes + (start_block)*block_size, 'bof');
         if status ~= 0
            error ('*** The value of starting block might be too large ***'); 
         end
      end
      
      nn = 0;		   % current plotting block data number
      nrem = 0;		% which block 10
      done = 0;		% not finish yet
      read_success = 1;
  else 
      % we are plotting real time and moving on to the next 50 blocks
      start_block = start_block + plot_blocks;           
   end % ~real_time && done
   if start_block-old_endblock<plot_blocks
       continued_flag=1; 
       skip_dist = (start_block-old_endblock)*pts_per_block;
   else
       skip_dist=0;
   end
end % while ~done

function fileName = getDataFileName;
% fileName = getDataFileName returns the file name of the data file.
% Preforms error checking.
%
% Fab, March 1998.

fileName='';
while isempty(fileName);   
   test_string = get_latest_file;
   fileName = input(['Enter data file name (default: ' test_string ' , ''q'' to quit): '], 's');
%   fileName = input('Enter data file name (default: ''testdata'', ''q'' to quit): ', 's');
   if strcmp(fileName, 'q')
      fclose('all');return;
 %     error(sprintf('\n****** File Loading aborted by user. Exiting ... *********'));
   elseif isempty(fileName)
      fileName = test_string;
   end
   if ~exist(fileName, 'file');
      fprintf('Can''t open file %s! Try again.\n', fileName);
      fileName = '';
   end
end

function fid = fileOpen(fileName);
% fid = fileOpen(fileName) returns the file ID for the file fileName.
% Preforms error checking.
%
% Fab, March 1998.

[fid, error_message] = fopen_odas(fileName, 'r');
if ~isempty(error_message), disp(error_message); end
if fid == -1,
   error(sprintf('Error opening file %s !\n', fileName));
end


function [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid)
% [nRow, nCol, nSlowCol, slowCh, fastCh] = load_ch_setup(fid) 
% loads the channel matrix from the file fileName.
% Returns the number of rows in the matrix, the number of columns in the matrix,
% the number of slow columns, a vector containing the slow channels,
% and a vector containing the fast channels.
% 
% NOTE: Consider making a standalone m file because this same function
%       gets used by plot_VMP.m
% 
% Replaces load_ch_setup.m written by L. Zhang.
% Fab, March 1998.
% AWS - 2010-01-14 changes for odas v6 and up

fseek(fid,0,'bof'); % Rewind to the beginning of the file
header = fread(fid, 64, 'ushort');%read the header
header_version = bitshift(header(11), -8) + bitand(header(11), 255) /1000;

if header_version >= 6
    nFast = header(29);
    nSlowCol = header(30);
    nRow  = header(31);
    nCol = nSlowCol + nFast;
    matrixSize = nRow * nCol; %use this to check what we get out of the setup file string
    setupfilestr = char(fread(fid, header(12), 'char'));
    
    if isempty(setupfilestr)
        error('failed to extract setup file string from first data record');
    end
    setupfilestr = setupfilestr';
    
    cfg = setupstr(setupfilestr);
    rows = setupstr(cfg, 'matrix', 'row[0-9]+');
    matrix = [];
    for row = rows
      values = textscan(row{1}, '%d16');
      matrix = vertcat(matrix, values{1}');
    end
    
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast
        error('Error building channel matrix: number of fast channels does not aggree.');
    end

else

    nFast = header(29);% extract info
    nSlowCol = header(30);
    nRow = header(31);

    % rebuild the channel matrix
    nCol = nSlowCol + nFast;
    matrixSize = nRow * (nSlowCol+nFast);
    matrix = fread(fid, matrixSize, 'ushort');
    matrix = reshape(matrix, nCol, matrixSize/nCol)';


    % build the slowCh vector
    slowCh = matrix(:, 1:nSlowCol)';
    slowCh = slowCh(:);
    if length(slowCh) ~= nSlowCol*nRow,
        error('Error building channel matrix: number of slow channels does not agree.');
    end

    %build the fastCh vector
    fastCh = matrix(1, nSlowCol+1:nCol)';
    if length(fastCh) ~= nFast,
       error('Error building channel matrix: number of fast channels does not aggree.');
    end
end
return


function [status, S, H] = read_block10(fid, S, H, whichbg, paramv)
% READ_BLOCK10
% support real-time plotting of data from a horizontal profiler using plot_HMP
%
% <latex>
% \index{Type A!setupstr}
%</latex>
%
% Function that supports the real-time plotting of data from a horizontal 
% profiler using plot_HMP.  It is not intended for direct use.  Reads 10 blocks
% of data at a time.
%
% |[status, S, H] = READ_BLOCK10(fid, S, H, whichbg, paramv)|
%
% * [in] _fid_ - open file ID with data to read
% * [in] _S_ - ?
% * [in] _H_ - header record from last time
% * [in] _whichbg_ - block group number (0 = first 10 blocks, 1 = second 10
%          blocks, etc.) max=4
% * [in] _paramv_ - vector of  input file specific parameters as below:
% 
% * [out] _status_ - 0 for success, -1 for failure
% * [out] _S_ - ?
% * [out] _H_ - updated header (if whichbg = 0)
%
% since more data gets read in than can be seen on a plot, the data for
% each channel is decimated with the minimum and maximum values from
% a user-specified time period replacing all of the individual readings
% for that time period. (this is specified in terms of points wanted
% per block)
%
% note that if all nblocks are not read successfully this returns a status
% flag indicating that it failed
%
% *Version History:*
% 
% * original version by Limin Zhang
% * 1995-02-01 (RAMO) modified extensively
% * 2004-05-07 (IG) modified to return a header on every block
% * 2012-04-30 (WID) Matlab plublishing support added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all "Magic" numbers here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size of integers from data file in bytes
size_of_integer=2;

% number of blocks to try and read when called
nblocks = 10;

% maximum number of successive failed reads allowed before we assume the input
% file isn't getting any bigger and its time to pack up (EOF)
max_failed_reads = 3;

% how long to pause before trying another read (in seconds)
onetick=2;

% paramv indices
header_size_i=1;        % index of the header size parameter 
rows_i =2;              %          number of rows of channels
cols_i =3;              %          number of columns of channels
num_slow_col_i=4;       %          number of slow channel columns
recs_per_block_i=5;     %          records per block
points_wanted_i=6;      %          points wanted (ie decimation factor)

% no magic numbers should appear below this line

% header size (number of integers in header)
header_size = paramv(header_size_i);

% header size (in bytes)
hsb = header_size*size_of_integer;

% rows of channels
rows = paramv(rows_i);

% columns of channels
cols = paramv(cols_i);

% number of columns of slow channels
num_slow_cols = paramv(num_slow_col_i);

% number of slow channels
num_slow_ch=num_slow_cols*rows;

% number of fast channels
num_fast_ch=cols-num_slow_cols;

% number of records per block
% one record is an array of (size rows x cols) slow/fast channel readings
recs_per_block=paramv(recs_per_block_i);

% total (for all channels) data per block (integers), excluding header
data_per_block = rows*cols*recs_per_block;

% total (for all channels) data per block (integers), including header
block_size = data_per_block+header_size;

% pointer to start of data within a block
sod = header_size+1;

% space required (integers) for nblocks blocks of slow channel data
slow_ch_len = num_slow_ch * recs_per_block * nblocks;

% space required for nblocks blocks of fast channel data
rows_per_block = rows*recs_per_block;

% number of readings for each slow channel per block
vals_per_slow_per_block = recs_per_block*nblocks;

%  points wanted per block for plotting (ie decimation factor)
%  one will be min the other max of a range of readings
points_wanted = paramv(points_wanted_i);

% total points wanted per channel
pts_per_ch = points_wanted*nblocks;

% pre-allocate storage needed

% for one block
blockin = zeros(1,block_size);

% for all the data
B = zeros(1, data_per_block*nblocks)+NaN;

% for the decimated data
T = zeros(1,pts_per_ch);
TT = zeros(1,pts_per_ch);
% reshape for min and max
TT = reshape(TT,2,pts_per_ch/2);

% blocks read in so far
blocks_read=0;

% Where the next block of data to be put in B starts
sob = 1;
% Where the next block of data to be put in B ends
eob = data_per_block;


% if input file is really at end of file
% (rather than stil being written as we are reading it)
EOF=0;

% if we are ready to read in another block and have not reached EOF
waiting=1;

% how many times in a row has a read failed to return a block of data
failed_reads = 0;


while waiting           % MAIN LOOP

% try to read a block;
  [blockin,count] = fread(fid, block_size,'short');
   
% did we get a full block?
  if (count==block_size)
% we got one
    failed_reads = 0;
    blocks_read=blocks_read+1;
% do we update the header ?
    if (blocks_read ==1) %& (whichbg == 0)
      H = blockin(1:header_size);
    end

% stuff the data into the big vector    
    B(sob:eob)=blockin(sod:block_size);    
% do we have enough?
    if (blocks_read >= nblocks)
      waiting = 0;
    else
% update pointers into B for the next incoming block
      sob=sob+data_per_block;
      eob=eob+data_per_block;
    end

  else
% if we didn't get a full block on the last read attempt
    failed_reads = failed_reads+1;
% do we call it quits?
    if (failed_reads > max_failed_reads)
      EOF=1;
      waiting=0;
    else
% or do we move back to where we were on the file and try again later
      fseek(fid, -count*size_of_integer, 'cof');
%     disp(['[',num2str(blocks_read),':',num2str(failed_reads),'] Waiting ',...
%         num2str(onetick),' seconds']);
      pause(onetick);
    end
  end
end   % of MAIN (reading loop) (because we got all nblocks blocks or hit EOF)

if EOF | (blocks_read < nblocks)
%signal failure
  status= -1;
else
%signal success
  status = 0;
  B = reshape(B,cols,rows_per_block*nblocks);
  SS = B(1:num_slow_cols,:);
  SS = SS(:);

  SS = reshape(SS,num_slow_ch,vals_per_slow_per_block);
  pts2 = recs_per_block/points_wanted; 
  if (pts2 ~=fix(pts2)) % this test was added by Fab, April 1998
     error('Sampling rate must be a power of 2.');
  end
  col = points_wanted * nblocks;
  aa = whichbg * pts_per_ch + 1;
  bb = (whichbg + 1) * pts_per_ch;
    for i=1:num_slow_ch
      T = SS(i,:)'';
      T = reshape(T, pts2*2, col/2);
      TD(1,:) = min(T);
      TD(2,:) = max(T);
      S(aa:bb,i) = TD(:);
    end
  pts2 = pts2 * rows;
[b,a]=butter(8,25/512);
  for i=1:num_fast_ch
        if (i>600);% Temporarily disabled by RGL, used to be (i>6)
          T = B(i+num_slow_cols,:);
          z=filtic(b,a,mean(T(1:32))*ones(1,8),mean(T(1:32))*ones(1,8));
          T=filter(b,a,T,z);
          B(i+num_slow_cols,:) = T;
        end
    T = reshape(B(i+num_slow_cols,:),pts2*2,col/2);
    TT(1,:) = min(T);
    TT(2,:) = max(T);
    S(aa:bb,i+num_slow_ch) = TT(:);
  end
end  % if success


function [tstr,pstr,hstr,sstr1,sstr2,day]=calibrations_plotting(S,ch_pars,mean_seq,mean_names,pts_per_block,day_flag,two_ch_var)
% ==========================================================================
%> @ingroup typeA
%> @file calibrations_plotting.m
%> @brief Internal function used by plot_HMP; Not to be called directly.
% ==========================================================================
%> @brief Sub-routine servicing plot_tomi.m, via plot_ch_data.m. Produces
%> calibrated average values of temperature (Sea-Bird), pressure, heading,
%> and speed, for inclusion on plot_tomi/plot_ch_data plots. Produces the
%> average values for the first block of whatever data is passed to the
%> routine.
%>
%> [tstr,pstr,hstr,sstr1,sstr2,day]=calibrations_plotting(S,ch_pars,mean_seq,mean_names,pts_per_block,day_flag,two_ch_var);
%>
%> @param S - data matrix
%> @param ch_pars - calibration coefficients, used to produce calibrated values
%>   for some of the variables listed in mean_seq/mean_names
%> @param mean_seq,mean_names - channel numbers, matrix positions, and channel
%>   names for the values to be printed on the figure (see plot_tomi.m)
%> @param pts_per_block - points per block
%> @param day_flag - determines whether or not the time of printing is shown
%> @param two_ch_var (optional) -  list of channels existing as even/odd pairs
%>   for a single variable (used for the new scount board).  If no input
%>   is provided, assume have only single-channel data.
%>
%> @retval Strings for temperature, pressure, heading, speed (two
%>   flowmeters) and day
%>
%> <b>Version History:</b>
%> - 2004-05-01 (IG) from code in plot_tomi.m (L. Zhang and others)
%>   Modification: Changed the flags used to locate the desired channels,
%>   now based on the names used in the main ODAS setup file
%>   Added code to deal with the new scount board; if the new board is in
%>   use (as reflected by the two_ch_var flag), the calibrations are done
%>   without wrap-around corrections, and without referring to additional
%>   frequency range parameters and frequency multiplication factors that
%>   were required with the old board.
%> - 2007-09-19 (RGL) added recognition of AE EM current meter and PNI
%>   compass.

if nargin<7, two_ch_var=[]; end

% Parameters
textFormat = '%7.2f'; % formatting string for average output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Locate the channels for the various output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the channels for which we must compute the
% average of the first block 
pres_ch = strmatch('Pres',mean_names);
temp_ch = strmatch('SBT2',mean_names);
u1_ch = strmatch('U',mean_names);
u2_ch = strmatch('V',mean_names);
ref_ch = strmatch('Cref',mean_names);
cos_ch = strmatch('Ccos',mean_names);
sin_ch = strmatch('Csin',mean_names);
mx_ch  = strmatch('Mx',mean_names);
my_ch  = strmatch('My',mean_names);

% Check which channels are associated with the new scount board (affects
% how they are calibrated below)
% IG, June 2004
if ~isempty(two_ch_var)
    if ~isempty(pres_ch)&&~isempty(find(two_ch_var(:,1)==mean_seq(pres_ch,1))), new_pres = 1; else new_pres=0; end
    if ~isempty(temp_ch)&&~isempty(find(two_ch_var(:,1)==mean_seq(temp_ch,1))), new_temp = 1; else new_temp=0; end
    if ~isempty(u1_ch)&&~isempty(find(two_ch_var(:,1)==mean_seq(u1_ch,1))), new_u1 = 1; else new_u1=0; end
    if ~isempty(u2_ch)&&~isempty(find(two_ch_var(:,1)==mean_seq(u2_ch,1))), new_u2 = 1; else new_u2=0; end
    if ~isempty(ref_ch)&&~isempty(find(two_ch_var(:,1)==mean_seq(ref_ch,1))), new_head = 1; else new_head=0; end
else
    new_pres=0; new_temp=0; new_u1=0; new_u2=0; new_head=0;
end

% Determine whether or not we have valid heading data with old FG compass
have_heading_FG = 0;
if (~isempty(cos_ch) && ~isempty(sin_ch) && ~isempty(ref_ch))
   have_heading_FG = 1;
end

% Determine whether or not we have valid heading data with new PNI
% magnetometer
have_heading_PNI = 0;
if (~isempty(mx_ch) && ~isempty(my_ch))
   have_heading_PNI = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the strings for the output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output average temperature if this was requested
if ~isempty(temp_ch) 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Formula for Sea-Bird probes
	%
	%   Temp = 1/{a + b[ln(f0/f] + c[sqr(ln(f0/f))] 
	%		+ d[cub(ln(f0/f)]} - 273.15
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IG, Jun. 2004 - Note: have changed the locations of fref and period,
    % to reflect the format in the main ODAS data file
    coef = ch_pars(temp_ch,:); % all coefficients for SeaBird Thermometer
	fref   = coef(6); % reference frequency on frequency-to-number converter of ODAS, usually 10e6Hz;
    period = coef(7); % p_count_set on ODAS board 4 to 256 in factors of 2;
    cons = fref * period;
    T = S(1:pts_per_block, mean_seq(temp_ch,2));
    T = T';
    if ~new_temp, T = T + 2^16; end % assume sign wrap-around on SeaBird counter
    T = cons ./ T;% SeaBird Thermometer frequency in Hz
    if ~new_temp
		w = find(T < 4000.0);% Check for wrap-around on ODAS f-to-num converter
		if ~isempty(w); % remove wrap-around by recalculating
			T(w) = S(w,mean_seq(temp_ch,2));
			T(w) = cons ./ T(w);
		end
    end
	T = log(coef(5) ./ T); 	% compute ln(f0/f)
	pcoef = fliplr(coef(1:4)); % reverse order of g h i j coefficients for SeaBird Thermometer 
	T = polyval(pcoef, T);
	T = 1 ./ T - 273.15;
	temp = sprintf(textFormat,mean(T)); 
else      
	% no temp ch given, leave it empty
	temp = blanks(8);
end
tstr = ['Temp = ', temp, ' ^oC'];

% output average pressure if this channel is assigned
if ~isempty(pres_ch)
	P = S(1:pts_per_block, mean_seq(pres_ch,2)); 
	coef = ch_pars(pres_ch,:);% Pressure coefficients
	coef = fliplr(coef); % reverse order for polyval
	P = polyval(coef, P); % Polynomial evaluation
	pres = sprintf(textFormat, mean(P));
else
	% not given, leave it empty
	pres = blanks(8);
end
pstr = ['Pres = ', pres, ' dBar'];

% compute and output heading direction
if have_heading_FG
%     MN = 21; 		% define true north
	HR = mean(S(1:pts_per_block, mean_seq(ref_ch,2))); %RGL changed this to compute mean first
	HS = mean(S(1:pts_per_block, mean_seq(sin_ch,2))); %and avoid wrap around problems
	HC = mean(S(1:pts_per_block, mean_seq(cos_ch,2))); 
	DEG = atan2(HS - HR, HC - HR) * 180 / pi;
	comp = sprintf(textFormat,DEG);
	hstr = ['Hding = ', comp, '^o'];
elseif have_heading_PNI
	My = mean(S(1:pts_per_block, mean_seq(my_ch,2))); %and avoid wrap around problems
	Mx = mean(S(1:pts_per_block, mean_seq(mx_ch,2))); 
	coef = ch_pars(mx_ch,:);% Mx coefficients
    Mx = (Mx - coef(1)) ./ coef(2); % In micro-tesla
	coef = ch_pars(my_ch,:);% My coefficients
    My = (My - coef(1)) ./ coef(2); % In micro-tesla
	DEG = atan2(My, Mx) * 180 / pi;
	comp = sprintf(textFormat,DEG);
	hstr = ['Hding = ', comp, '^o'];
    
else
  % not wanted, set to empty string
  hstr = ['Hding = ', blanks(8), '^o'];
end

% This is the old calculation for a propellor current meter, which is no
% longer supported 
% compute speed U1 and output it 
% if ~isempty(u1_ch) 
%      U1 = S(1:pts_per_block, mean_seq(u1_ch,2)); 
%      coef = ch_pars(u1_ch,:); % all coefficients for propelor anemomter
%      fref   = coef(7); % reference frequency on frequency-to-number converter of ODAS, usually 10e6Hz;
%      period = coef(8); % p_count_set on ODAS board 4 to 256 in factors of 2;
%      scale  = coef(6); % Divisor in phase-locked-loop of propelor electronics (only for old scount board)
%      U1 = U1'; % make into row vector
%      if ~new_u1
%          w = find (U1 <=0);% Check for wrap around on ODAS freq-to-num converter
%          if ~isempty(w); U1(w) = U1(w) + 2^16; end % add offset for wrap-around
%      end
%      cons = fref * period;
%      U1 = cons ./ U1;% Propelor frequency in Hz
%      if ~new_u1, U1 = U1/scale; end% Adjust for PLL in electronics
%      pcoef = fliplr(coef(1:4)); % reverse order of calibration coefficients for polyval 
%      U1 = polyval(pcoef, U1);
%      U1 = U1(find(abs(U1 - mean(U1)) < 3*std(U1))); % remove outliers
%      if ~isempty(U1);
%         U1 = mean (U1);
%      else
%         U1 = -99;
%      end
%      u1str = sprintf(textFormat, U1);
%      sstr1 = ['U1 = ', u1str, ' m/s'];
% else
%       % not wanted, set to empty string
%       sstr1 = ['U1 = ',blanks(8), ' m/s'];
% end 
%
% This is the new calculation for an Alec Electronics Ellectro-magnetic
% current meter
% compute speed U and output it 
if ~isempty(u1_ch) 
	U = mean(S(1:pts_per_block, mean_seq(u1_ch,2)));
    coef = ch_pars(u1_ch,:); % all coefficients for EMCM
	coef = fliplr(coef); % reverse order for polyval
	U = polyval(coef, U); % Polynomial evaluation
    u1str = sprintf(textFormat, U);
    sstr1 = ['U = ', u1str, ' m/s'];
else
      % not wanted, set to empty string
      sstr1 = ['U = ',blanks(8), ' m/s'];
end 

% This is the new calculation for an Alec Electronics Ellectro-magnetic
% current meter
% compute speed V and output it 
if ~isempty(u2_ch) 
	V = mean(S(1:pts_per_block, mean_seq(u2_ch,2)));
    coef = ch_pars(u2_ch,:); % all coefficients for EMCM
	coef = fliplr(coef); % reverse order for polyval
	V = polyval(coef, V); % Polynomial evaluation
    u2str = sprintf(textFormat, V);
    sstr2 = ['V = ', u2str, ' m/s'];
else
      % not wanted, set to empty string
      sstr2 = ['V = ',blanks(8), ' m/s'];
end 


% output the processing (current) date
% only at 5th ten block segment or all done
if day_flag
      % get current date
      d = date;

      % get current time
      clk = clock;
      t = sprintf('%2.0f:%2.0f:%2.0f', clk(4),clk(5),clk(6));
      if t(7)==' '
         t(7)='0';
      end

      % put all of them together and output date string
      day = ['Printed on: ',d, ', ',t];
  else 
      day = [blanks(8)];
end


function [plot_seq] = assign_plot_seq(ch_plot_seq, slow_ch, fast_ch)
%> @defgroup typeA Type A functions
% ===========================================================================
%> @ingroup typeA
%> @file assign_plot_seq.m
%> @brief Internal function used by plot_HMP; Not to be called directly.
% ===========================================================================
%> @brief Assign each of 16 plotting channels the corresponding location
%> in the matrix. This matrix location index enable us to  
%> extract a particular channel data from records of blocks.
%>
%> Records is the smallest repeated unit in the block, i.e.,
%> one matrix of data is treated as a record of data
%>
%> Channel Matrix Location Index:
%> \code
%> 1 (k+1) (k+2) (k+3) (k+4) ... (k+j)
%> 2 (k+1) (k+2) (k+3) (k+4) ... (k+j)
%> 3 (k+1) (k+2) (k+3) (k+4) ... (k+j)
%>     .  .  . 
%>     .  .  . 
%>
%> k (k+1) (k+2) (k+3) (k+4) ... (k+j)
%> \endcode
%>
%> Here, 'k' is the number of slow channels and 'j' is the 
%> number of fast channels.
%>
%> function [plot_seq] = assign_plot_seq(ch_plot_seq, slow_ch, fast_ch)
%> 
%> @param ch_plot_seq
%> @param slow_ch
%> @param fast_ch
%>
%> @retval plot_seq - a matrix of two columns; The first
%> column contains the channel numbers to be plotted and 
%> the second column is the channel matrix location index. 

% the maxinum number of processing channels specified in 'ch_data.txt'
max_ch = length(ch_plot_seq);

% preallocat storage
seq = (1:1:max_ch);

% the number of slow channels
k = length(slow_ch);		

% assign the channel matrix location to each of these channels 
for i=1:max_ch
  % is ch_plot_seq number valid?
  if (any(slow_ch == ch_plot_seq(i)) | any(fast_ch == ch_plot_seq(i))) 
    j = find(slow_ch == ch_plot_seq(i));	% is it slow channel?

    if ~isempty(j), seq(i) = j(1); end		% first matrix location index,

    j = find(fast_ch == ch_plot_seq(i));	% is it fast channel?

    if ~isempty(j), seq(i) = k + j(1); end	% first matrix location index
    					
  else
     msg = ['*** Error: invalid channel number: ', num2str(ch_plot_seq(i))];
     error(msg);
  end
end

% construct two columns: 
% 1st = channel number; 2nd = matrix location index
plot_seq = [(ch_plot_seq)'; seq]';
return


function plot_ch_data(fileName,S,H,nn,nrem,seq,ch_names,mean_seq,mean_names,ch_pars,pts_per_block,pts_per_ch,start_block,plot_blocks,day_flag,two_ch_var)
%==========================================================================
%> @ingroup typeA
%> @file plot_ch_data.m
%> @brief Internal function used by plot_HMP.m; do not call directly.
%==========================================================================
%> @brief Internal function used by plot_HMP.m (formerly called
%>  plot_tomi.m). Given various input variables, plot the data, print out 
%>  some average values.
%>
%> plot_ch_data(fileName,S,H,nn,nrem,seq,ch_names,mean_seq,mean_names,ch_pars,
%>  pts_per_block,pts_per_ch,start_block,plot_blocks,day_flag,two_ch_var)
%>
%> @param fileName - name of the file
%> @param S - data matrix
%> @param H - header data
%> @param nn, nrem - number of reads and number (0-4) of 10-block groups read 
%>  thusfar; defined in plot_HMP.m
%> @param seq, ch_names - channel numbers, matrix positions, and channel names 
%>  for the line plots (see plot_HMP.m)
%> @param mean_seq, mean_names - channel numbers, matrix positions, and channel
%>  names for the values to be printed on the figure (see plot_HMP.m)
%> @param ch_pars - calibration coefficients, used to produce calibrated values
%>  for some of the variables listed in mean_seq or mean_names
%> @param pts_per_block, pts_per_ch - points per block and channel
%> @param start_block,plot_blocks - starting block and number of blocks to be
%>  plotted (used to determine the x-axis)
%> @param day_flag - determines whether or not the time of printing is shown
%> @param two_ch_var - (optional) list of channels existing as even and odd
%>  pairs for a single varialbe (used for the new scount board).  If no input
%>  is provided, assume have only single-channel data.
%> 
%> Routines called: calibrations_plotting, sub-routines also included below.
%>
%> <b>Version History:</b>
%> - from plot_tomi.m (by L. Zhang and others)  
%> - 2004-06-01 (IG) Modifications: Changed the figure sizing slightly to
%>  fit on the screen. Screen now updated every 10 blocks, 50-block 'window'
%>  slides along the data so that the most recent data is always at the right
%>  edge. Calibrated values, shown as text, now represent the values in the most
%>  recent block, rather than the earliest one shown on the plot (as in the
%>  original version of plot_tomi). Added a section of code to deal with the
%>  new scount board, in which data must be combined from odd and even channels.
%> - 2005-01-01 (IG) modified handling of split channels

if nargin<16, two_ch_var=[]; end

% Size of the figure on the screen and on the paper (normalized units on the screen, inches on paper)
figureHeight = 0.875;                       % Position so fill the screen, same aspect ratio as on paper
figureWidth = (8/10.5)*figureHeight;    
figurePos = [0.1 0.075 figureWidth figureHeight];
figureUnits = 'normalized';
figureTag = 'plot_tomi';
paperPos = [0.25 0.25 8 10.5];
paperUnits = 'inches';

% Other plotting parameters
inc = 1/20;% global parameters for plotting 
axisHeight = 1/26;
axisWidth = 0.80;
axisXPos = 0.12;
topLineYPos = .985; % coordinates for the average channel readings (w.r.t. hTextAxis)
aveLineYPos1 = 0.53;
aveLineYPos2 = 0.49;
lineTab1 = 0.1;
lineTab2 = 0.35;
lineTab3 = 0.6;

% initialize names of months
Months = ['Jan'; 'Feb'; 'Mar'; 'Apr';
	       'May'; 'Jun'; 'Jul'; 'Aug';
	       'Sep'; 'Oct'; 'Nov'; 'Dec'];
       
% If we are dealing with the new scount board (i.e., have two channels for
% some variables), then find and combine the appropriate channels, modify
% the lists of channels to plot etc. accordingly
% If not using the new scount board, this will do nothing.
% Added by IG, June 2004
max_ch = size(seq,1);       % Number of channels to plot
if ~isempty(two_ch_var)     % If some channels are odd/even parts of single data variables, combine as necessary
    for ii=1:size(two_ch_var,1)
        junk_ii = [find(seq(:,1)==two_ch_var(ii,1)) find(seq(:,1)==two_ch_var(ii,2))];
        if length(junk_ii)==2
            matrix_pos = seq(junk_ii,2);
            max_ch=max_ch-1;
            seq(junk_ii(2),:)=[];
            junk = deblank(ch_names(junk_ii(1),:)); ch_names(junk_ii(1),1:length(junk))=[junk(1:end-1) ' '];
            ch_names(junk_ii(2),:)=[];
            junk1 =find(mean_seq(:,1)==two_ch_var(ii,1)); junk2 = find(mean_seq(:,1)==two_ch_var(ii,2));
            if ~isempty(junk1)
                mean_seq(junk2,:)=[];
                junk = deblank(mean_names(junk1,:)); mean_names(junk1,1:length(junk))=[junk(1:end-1) ' '];
                mean_names(junk2,:)=[];
                ch_pars(junk2,:)=[];
            end
        else
            junk_ii = [find(mean_seq(:,1)==two_ch_var(ii,1)) find(mean_seq(:,1)==two_ch_var(ii,2))]; 
            matrix_pos = mean_seq(junk_ii,2);
            mean_seq(junk_ii(2),:)=[];
            junk = deblank(mean_names(junk_ii(1),:)); mean_names(junk_ii(1),1:length(junk))=[junk(1:end-1) ' '];
            mean_names(junk_ii(2),:)=[];
            ch_pars(junk_ii(2),:)=[];
        end
        junk = find(S(:,matrix_pos(1))<0); S(junk,matrix_pos(1))=S(junk,matrix_pos(1))+2^16;
        junk = find(S(:,matrix_pos(2))<0); S(junk,matrix_pos(2))=S(junk,matrix_pos(2))+2^16;
        S(:,matrix_pos(1)) = S(:,matrix_pos(1))+S(:,matrix_pos(2))*2^16;
        S(:,matrix_pos(2)) = NaN;       % As a precaution, remove the second part of the values to prevent them being used by accident
    end     
end

% If we are just starting, close all open figures and generate one of the
% correct size; otherwise just clear the figure.
if isempty(findobj('type','figure','tag',figureTag))
    close all;
    hFigure = figure('units',figureUnits,'position',figurePos,'paperunits',paperUnits,'paperposition',paperPos,'tag',figureTag);
else
    clf; % clear previous figures
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Plot Each Channel Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting from bottom to top
for i=max_ch:-1:1     
	if i <= 8; % set y position of upper left corner of the axis
		axisYPos = 0.12 + (17 - i)*inc;
	else 
		axisYPos = 0.02 + (17 - i) * inc;
	end

	SV = S(:,seq(i,2));% SV contains a copy of plotting channel data
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% NOTE: computing 'ss' is necessary, otherwise NaN is
	%       always returned as the result of the following
	%       'min' and 'max' function calls.
	%
	% I do not believe that this is true, RGL
    % IG - does not in fact appear to be unnecessary
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ss = (nrem+1)*pts_per_block*10;
	ymin = min(SV(1:ss));
	ymax = max(SV(1:ss));
 
 
	% create an axis, data plotting occurs in this axis
	hPlotAxis = axes('position',[axisXPos axisYPos axisWidth axisHeight]);
	
	% create abscissa vector
	abscissa = linspace(start_block, start_block+plot_blocks, pts_per_ch);
	
	plot(abscissa, SV); % plot the current vector
	
	if ymin == ymax; % check if ymin equals ymax
		ymin = ymin - 1;
		ymax = ymax + 1;
	end

	d = (ymax - ymin)*0.10; % this is for output 'y' min and max labels in right position
	
	% set the range and ticks as above
	set(hPlotAxis, 'ylim', [ymin ymax],'YTickLabel', [],'Ytick',[]); % RGL added 'Ytick'
	ymaxString = sprintf('%.0f ', ymax);
	yminString = sprintf('%.0f ', ymin);
	hYmaxString = text(start_block, ymax-d, ymaxString, 'HorizontalAlignment', 'right','fontsize',10);
	hYminString = text(start_block, ymin+d, yminString, 'HorizontalAlignment', 'right');
	
	% set the range of x values and tick mark locations
	set(hPlotAxis, 'xlim', [abscissa(1) abscissa(pts_per_ch)],...
		'XTick', [start_block:(start_block+plot_blocks)], ...
		'XTickLabel', xTickLabelString(start_block, plot_blocks));

	% only two axes with x labels; 8th and last from the bottom 
	if i ~= 8 & i ~= max_ch 
		set(hPlotAxis, 'XTickLabel' , []);
	end
	if i == max_ch,
		xlabel('block number')
	end

	% set other axis properties
	set(hPlotAxis, 'FontName', 'Helvetica', ...
		'FontSize', 10,...
		'ticklength', 0.5*get(hPlotAxis,'TickLength'));
	
	% output axis labels for y axis
	hTextString=text(start_block, (ymax+ymin)/2, ch_names(i,:),...
		'HorizontalAlignment', 'right',...
		'FontName', 'Helvetica',...
		'FontSize', 10);
end % for loop
  
% create an invisible axis for text output (all text is placed relative to this)
hTextAxis = axes('position', [0 0 1 1],...
	'visible', 'off');
  
% output data file name at the upper left corner
setTextParameters(text(lineTab1, topLineYPos, ['File: ' strrep(fileName,'_','\_')]));

% output measuring time at the centre of the plot
d = sprintf('%2d-%s-%4d, ',H(6),Months(H(5),:),H(4));
if d(1)==' ';
	d(1)='0';
end
t = sprintf('%2d:%2d:%2d UT', H(7), H(8), H(9));
mtime = ['Measured on ', d, t];
setTextParameters(text(lineTab3, aveLineYPos1, [mtime]));
  
% Output the desired mean, calibrated values
[tstr,pstr,hstr,sstr1,sstr2,day]=calibrations_plotting(S(end-pts_per_block+1:end,:),ch_pars,mean_seq,mean_names,pts_per_block,day_flag,two_ch_var);

setTextParameters(text(lineTab1, aveLineYPos1, [tstr]));
setTextParameters(text(lineTab3, aveLineYPos2, [pstr]));
setTextParameters(text(lineTab1, aveLineYPos2, [hstr]));
setTextParameters(text(lineTab2, aveLineYPos1, [sstr1]));
setTextParameters(text(lineTab2, aveLineYPos2, [sstr2]));
setTextParameters(text(lineTab3, topLineYPos, [day]));

% flush and draw now
drawnow;

function s=xTickLabelString(startBlock, numBlocks)
% string = xTickLabelString(block, plot_blocks) returns a cell array string of 
% abscissa labels spaced 5 blocks apart, starting at startBlock, ending at 
% startBlock+numBlocks.
% Fab, March 1998.

tickVector = ones(numBlocks+1,1)*NaN;
markers = startBlock:5:startBlock+numBlocks; % the values at the tick marks
index = markers-markers(1);
tickVector(index+1) = markers;
s = '';
for i = 1:numBlocks+1
   if ~isnan(tickVector(i))
      s = strvcat(s,num2str(tickVector(i)));
   else
      s = strvcat(s, blanks(1));
   end
end
return;

function setTextParameters(textHandel)
% setTextParameters(textHandle) sets the Fontname and Fontsize parameters
% for the output text string with handle textHandel
% Fab, March 1998.
set(textHandel, 'FontName', 'Helvetica', 'FontSize', 10);


function [ch_plot_seq,ch_names,ch_mean_seq,ch_mean_names,ch_mean_coef,two_ch_var]=load_setup_data(fname)
%===========================================================================
%> @ingroup typeA
%> @file load_setup_data.m
%> @brief Internal function used by plot_HMP; do not use directly
%===========================================================================
%> @brief A service function to support real-time plotting of data collected
%>  with a horizontal profiler using plot_HMP and  not intended for direct use.
%>
%>
%> Sub-routine servicing plot_tomi.m. Given the name of an ODAS setup file,
%> extracts the information required by plot_tomi.m for plotting etc.  Note
%> that the file depends on the use of particular keywords, which must be
%> present in the text file and defined as expected by the routine (defined
%> below).
%>
%> [ch_plot_seq,ch_names,ch_mean_seq,ch_mean_names,ch_mean_coef,two_ch_var]=
%>   load_setup_data(fname)
%>
%> @param fname - name of the setup file (a .txt file)
%>
%> @retval ch_plot_seq - channel numbers to be plotted by plot_tomi.m
%> @retval ch_names - names of these channels
%> @retval ch_mean_seq - channel numbers for which average values may be 
%>  computed
%> @retval ch_mean_names - names of these channels
%> @retval ch_mean_coef - calibration coefficients for these channels
%> @retval two_ch_var - data variables split up over two channels (gives both 
%>  the first and the second channel that must be combined). This variable is
%>  optional; if used as an output variable and there are no two-channel data,
%>  returns an empty vector
%> 
%> Keywords required in the setup file:
%>  - channel - \<channel number\>,\<channel name\>,\<series of calibration
%>    coefficients\>
%>  - plotting - series of channel numbers to plot
%>  - plotaverages - series of channel numbers for which the average may be
%>       computed (kept separate from plotting, all elements of which must
%>       be plotted
%>  - averagenames - same as 'plotting', but for the averaging channels 
%>    Elements within each line must be comma-delimited
%>   
%> <b>Version History:</b>
%> - 2004-06-01 (IG) initial

if nargin<1, fname=[]; end

% Parameters
num_coeffs = 8;         % Number of coefficients expected by plot_tomi; ignore any extra coefficients in the setup file

% Try to open the setup file (currently allows for one mistake in entering the name)
if ~isempty(fname), fid = fopen(fname);end    
if isempty(fname) || fid==-1
    fname = input('Enter the setup file name: ','s');
    fid = fopen(fname);
    if fid==-1, error('Could not find the specified file'); end
end

% Extract the required values from the file.  Uses keywords to locate the
% desired lines, ignoring lines that are commented out with a hash-mark (#)
done = 0;
ch_num=[];
ii=1;
while ~done
    junk=fgetl(fid);
    if junk==-1, done=1; 
    elseif (length(junk)>1&&~strcmp(junk(1),'#'))        % ignore commented lines
        start_ii = findstr(junk,':')+1;
        delim_ii = [findstr(junk,',') length(junk)+1];
        if ~isempty(findstr(junk,'channel'))            % get the calibration data
            try ch_num(ii) = str2num(deblank(junk(start_ii:delim_ii(1)-1)));
            catch
                disp(''); 
            end
            for jj=1:length(delim_ii)-2
                ch_cal(ii,jj) = str2num(junk(delim_ii(jj+1)+1:delim_ii(jj+2)-1)); 
            end
            ii=ii+1; 
        elseif ~isempty(findstr(junk,'plotting'))       % get channel numbers to plot
            ch_plot_seq(1) = str2num(junk(start_ii+1:delim_ii(1)-1));
            for jj=1:length(delim_ii)-1
                ch_plot_seq(jj+1)=str2num(junk(delim_ii(jj)+1:delim_ii(jj+1)-1)); 
            end
        elseif ~isempty(findstr(junk,'plotaverages'))   % get channel numbers for averaging
            ch_mean_seq(1) = str2num(junk(start_ii+1:delim_ii(1)-1));
            for jj=1:length(delim_ii)-1
                ch_mean_seq(jj+1)=str2num(junk(delim_ii(jj)+1:delim_ii(jj+1)-1)); 
            end
        elseif ~isempty(findstr(junk,'plotnames'))      % get names of plotting variables
            start_ii = max(findstr(junk(start_ii+1:delim_ii(1)-1),' '))+start_ii;
            ch_names{1} = deblank(junk(start_ii+1:delim_ii(1)-1));
            for jj=1:length(delim_ii)-1
                ch_names{jj+1}=deblank(junk(delim_ii(jj)+1:delim_ii(jj+1)-1)); 
            end
        elseif ~isempty(findstr(junk,'averagenames'))   % get names of averaging variables
        start_ii = max(findstr(junk(start_ii+1:delim_ii(1)-1),' '))+start_ii;
        ch_mean_names{1} = deblank(junk(start_ii+1:delim_ii(1)-1));
        for jj=1:length(delim_ii)-1
            ch_mean_names{jj+1}=deblank(junk(delim_ii(jj)+1:delim_ii(jj+1)-1)); 
        end
        end
    end
end    
% Re-organize variables a bit, extract only the required calibration data
ch_plot_seq=ch_plot_seq';
ch_mean_seq = ch_mean_seq';
ch_mean_coef = zeros(length(ch_mean_seq),num_coeffs)*NaN;
for ii=1:length(ch_mean_seq)
    jj = find(ch_num==ch_mean_seq(ii));
    if ~isempty(jj), ch_mean_coef(ii,:) = ch_cal(jj,1:num_coeffs); end
end
clear ch_cal
ch_names = strvcat(ch_names);
ch_mean_names = strvcat(ch_mean_names);

% Deal with data that are spread over two channels (required for the
% new-scout board)
% Added by IG, June 2004
two_ch_var = [];
for ii=1:size(ch_names,1)-1
    junk1=deblank(ch_names(ii,:)); junk2=deblank(ch_names(ii+1,:));
    if strcmp(junk1(1:end-1),junk2(1:end-1))&&(strcmp(junk1(end),'E')||strcmp(junk2(end),'E'))&&(strcmp(junk2(end),'O')||strcmp(junk1(end),'O'))
        if strcmp(junk1(end),'E')
            two_ch_var = [two_ch_var; ch_plot_seq(ii) ch_plot_seq(ii+1)]; 
        elseif strcmp(junk2(end),'E')
            error(['Odd channel of ' deblank(ch_names(ii,:)) ' seems to be numbered lower than the even channel; there may be an error in the setup file']); 
        end
    end
end
for ii=1:size(ch_mean_names,1)-1
    junk1=deblank(ch_mean_names(ii,:)); junk2=deblank(ch_mean_names(ii+1,:));
    if strcmp(junk1(1:end-1),junk2(1:end-1))&&(strcmp(junk1(end),'E')||strcmp(junk2(end),'E'))&&(strcmp(junk2(end),'O')||strcmp(junk1(end),'O'))
        if ~isempty(two_ch_var)&&~isempty(find(two_ch_var(:,1)==ch_mean_seq(ii))), 
            return; end
        if strcmp(junk1(end),'E')
            two_ch_var = [two_ch_var; ch_mean_seq(ii) ch_mean_seq(ii+1)];
        elseif strcmp(junk2(end),'E')
            error(['Odd channel of ' deblank(ch_names(ii,:)) ' seems to be numbered lower than the even channel; there may be an error in the setup file']); 
        end
    end
end

% Close the file
fclose(fid);

if nargout<6, clear two_ch_var; end


function [start_block,end_block,nb,quit_flag]=input_blocks(real_time,done,plot_blocks,default_start);
%===========================================================================
%> @ingroup typeA
%> @file load_setup_data.m
%> @brief Internal function used by plot_HMP; do not use directly
%===========================================================================
%> @brief A service function to support real-time plotting of data collected
%>  with a horizontal profiler using plot_HMP and  not intended for direct use.
%>

%> Sub-routine servicing plot_tomi.m.  Asks the user to input start and/or
%> end blocks (depending on the type of plotting being done), checks the
%> numbers, and puts out value for use by plot_tomi.
%>
%>[start_block,end_block,nb,quit_flag]=input_blocks(real_time,done,plot_blocks,
%>   default_start);
%> 
%> @param real_time, done - flags indicating whether or not plotting is
%>     real-time, and whether or not we are done plotting a set of blocks
%> @param plot_blocks - number of blocks to plot in a window
%> @param default_start - default start_block value (usually 0 or previous
%>     end_block)
%>
%> @retval start_block, end_block - start and end blocks to be read from 
%>  TOMI file
%> @retval nb - number of sets of blocks; number of multiples of ten blocks
%>     (smallest number past end_block)
%> @retval quit_flag - set to 1 if the user has decided not to continue plotting
%>       (0 by default)
%>
%> <b>Version History:</b>
%>  - 2004-06-01 (IG) from code in plot_tomi.m (L. Zhang and others)

% Default values for the starting block and continuing plotting
if nargin<3, default_start=0; end
quit_flag=0;

% Determine the starting block.  Ask the user whether or not to continue
% plotting of we've already reached the end of a set of blocks, then ask
% for the starting block (if appropriate)
if ~real_time & done
    answer = input('Continue plotting (n/y): ','s');
    if strcmp(answer,'y')
        msg = ['Enter starting blocks (ENTER for default: ', int2str(default_start),'): '];
        start_block = input(msg);
    else 
        quit_flag=1;
    end
else
    start_block = input('Enter the number of blocks to skip (default: 0): ');
end

% If we are to continue plotting (most cases), check the starting block
% entered above, then ask for the ending block if appropriate (i.e., if not
% real time)
if ~quit_flag
	if isempty(start_block)% no input, set default to the beginning of the file
        start_block = default_start;
	elseif start_block < 0
       start_block = 0;
	end
	start_block = start_block - rem(start_block, 5);% Start on 5-block boundary
	if real_time == 0 % not real time plotting, then get the ending block
      end_block = start_block + plot_blocks;
      msg = ['Enter ending block (ENTER for default: ', int2str(end_block),'): ']; 
      end_block = input(msg);
	
      if isempty(end_block); % no input or invalid input, use default
         end_block = start_block + plot_blocks;
      elseif end_block < start_block
         end_block = start_block + plot_blocks;
      end
	
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      % NOTE:
      %      nb is the number of multiples of 10 blocks, 
      %      which is the smallest number bigger than 
      %      the end_block.                                     
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      nb = ceil((end_block - start_block)/10);
	else
        end_block = []; nb=[];
	end
else start_block=[]; end_block=[]; nb=[];   % Default case, just to keep the routine going if have quit
end

if nargout<4, clear quit_flag; end
