function [UFILE] = ufile_read_nbi(filename, plotyn)
% Read and plot the UFILEs produced by MC3

% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29210/P29210.NBI';

% INPUT 

% filename	name of the UFILE to read
% plotyn	keyword = 1 to plot, 0 otherwise
% tslice	time slice for plot 1D profile from 2D, set to 0 for no slice
%
% OUTPUT
% UFILE		structure containing the relevant data
%
% Example:
% UFILE = ufile_read_nbi(filename = 'P29880.NBI',plotyn = 1);
% UFILE = ufile_read_nbi(filename = 'P29880.NBI',plotyn = 1);

UFILE.filename = filename;

fid = fopen (filename, 'r');
k = 0;
do
  k = k + 1;
  fileheader{k} = fgetl(fid);
until (k == 11)
fclose(fid);
UFILE.header = fileheader;


% Read the 1st line
txt = strsplit(fileheader{1});
UFILE.shot =  txt{1,2};
UFILE.device = txt{1,3};
UFILE.data_dimension = str2num(txt{1,4});
UFILE.datum = txt{1,10};
clear txt
tslice=0.252;


% Read the X data label
if (UFILE.data_dimension == 2)
  read_line_no = 6;
  UFILE.x_label = fileheader{read_line_no}(1:30);
endif


% Read the signal name
read_line_no = 8;
txt = strsplit(fileheader{read_line_no});
UFILE.signal = fileheader{read_line_no}(1:30);
clear txt

% Read the time dimension
read_line_no = 10;
txt = strsplit(fileheader{read_line_no});
UFILE.M = str2num(txt{1,2});
clear txt

% Read the data dimension
read_line_no = 11;
txt = strsplit(fileheader{read_line_no});
UFILE.N = str2num(txt{1,2});
clear txt



% Size of the data array
UFILE.data_size = UFILE.M * UFILE.N;


% Now read the data
fid = fopen (filename, 'r');
number_of_lines_to_skip = 11;
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

% Find the time slice
if (tslice != 0)
  UFILE.tslice_index = max(find(UFILE.t <= tslice));
  UFILE.tslice = UFILE.t(UFILE.tslice_index);
  UFILE.tslice_data =  UFILE.data(UFILE.tslice_index,:);
endif



% plot the data
if (plotyn == 1)
  if (UFILE.data_dimension == 1)
    figure(1)
    plot(UFILE.t, UFILE.data, 'k-', 'linewidth', 2)
    title([UFILE.shot ' ' UFILE.device ' ' UFILE.datum ' ' UFILE.signal])
    xlabel('time (s)')
    ylabel(UFILE.signal)
    axis tight
    as = axis;
    axis ([as(1) as(2) as(3) 1.1*as(4)])
    if (tslice != 0)
      text(0.1*as(2), 1.05*as(4), ['t = ' num2str(UFILE.tslice) ' s, Value = ' num2str(UFILE.tslice_data, '%e')]);	
    endif

  elseif (UFILE.data_dimension == 2)
    [t,x] = meshgrid(UFILE.t, UFILE.x);
    figure(2, 'position', [100 100 1000 400])
    subplot(1,2,1)
    pcolor(t, x, UFILE.data')
    shading flat
    xlabel('time (s)')
    ylabel(deblank(UFILE.x_label))
    colorbar
    axis tight
    title([UFILE.shot ' ' UFILE.device ' ' UFILE.datum ' ' UFILE.signal])


    subplot(1,2,2)
    surf(t, x, UFILE.data')
    shading flat
    xlabel('time (s)')
    ylabel(deblank(UFILE.x_label))
    zlabel(UFILE.signal)
    title([UFILE.shot ' ' UFILE.device ' ' UFILE.datum ' ' UFILE.signal])
    
     if (tslice != 0)
      figure(3)
      plot(UFILE.x, UFILE.tslice_data, 'k-', 'linewidth', 2)
      title([UFILE.shot ' ' UFILE.device ' ' UFILE.datum ' ' UFILE.signal])
      xlabel(deblank(UFILE.x_label))
      ylabel(UFILE.signal)
      axis tight
      as = axis;
      text(0.1*as(2), 0.95*as(4), ['t = ' num2str(UFILE.tslice) ' s']);	
    endif
    
  endif
endif
