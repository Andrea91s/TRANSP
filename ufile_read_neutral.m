function [UFILE] = ufile_read_neutral(filename, plotyn)
% Read and plot the UFILEs produced by MC3
% 
% INPUT 
% filename	name of the UFILE to read
% plotyn	keyword = 1 to plot, 0 otherwise

%
% OUTPUT
% UFILE		structure containing the relevant data
%
% Example:
% filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U01/F29904.RCY';
% [UFILE] = ufile_read_neutral(filename, plotyn = 1);


UFILE.filename = filename;

fid = fopen (filename, 'r');
k = 0;
do
  k = k + 1;
  fileheader{k} = fgetl(fid);
until (k == 7)
fclose(fid);
UFILE.header = fileheader;


% Read the 1st line
txt = strsplit(fileheader{1});
UFILE.shot =  txt{1,2};
UFILE.device = txt{1,3};
UFILE.data_dimension = str2num(txt{1,4});
UFILE.datum = txt{1,10};
clear txt
read_line_no = 4;


% Read the X data label
if (UFILE.data_dimension == 2)
  read_line_no = read_line_no + 1;
  UFILE.x_label = fileheader{read_line_no}(1:30);
endif


% Read the signal name
read_line_no = read_line_no + 1;
txt = strsplit(fileheader{read_line_no});
UFILE.signal = fileheader{read_line_no}(1:30);
clear txt

% Read the time dimension
read_line_no = read_line_no + 2;
txt = strsplit(fileheader{read_line_no});
UFILE.M = str2num(txt{1,2});
clear txt

% Read the data dimension
if (UFILE.data_dimension == 1)
  UFILE.N = 1;
elseif (UFILE.data_dimension == 2)
  read_line_no = read_line_no + 1;
  txt = strsplit(fileheader{read_line_no});
  UFILE.N = str2num(txt{1,2});
  clear txt
end


% Size of the data array
UFILE.data_size = UFILE.M * UFILE.N;


% Now read the data
fid = fopen (filename, 'r');
number_of_lines_to_skip = read_line_no;
fskipl(fid, number_of_lines_to_skip);
UFILE.t = fscanf(fid, '%f', UFILE.M);
if (UFILE.data_dimension == 2)
  UFILE.x = fscanf(fid, '%f', UFILE.N);
endif
UFILE.data = fscanf(fid, '%f', UFILE.M * UFILE.data_size);
UFILE.data = reshape(UFILE.data, UFILE.M, UFILE.N);
fclose(fid);

% Read the last line
fid = fopen (filename, 'r');
nl = fskipl(fid, Inf);  
fclose(fid);
fid = fopen (filename, 'r');
fskipl(fid, nl - 1);
UFILE.bottomer = fgetl(fid);
fclose(fid);


% plot the data
if (plotyn == 1)
  if (UFILE.data_dimension == 1)
    figure
    plot(UFILE.t, UFILE.data, 'k-', 'linewidth', 2)
    title([UFILE.shot ' ' UFILE.device ' ' UFILE.datum ' ' UFILE.signal])
    xlabel('time (s)')
    ylabel(UFILE.signal)
    axis tight
    as = axis;
    axis ([as(1) as(2) as(3) 1.1*as(4)])
  
  endif
endif
