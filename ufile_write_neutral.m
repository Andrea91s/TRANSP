function [UFILE] = ufile_write_neutral(UFILE, filename)
% Write a UFILE
% 
% INPUT 
% UFILE		structure containing the relevant data to be written
% filename	name of the file to be saved

% Example:
% UFILE = ufile_read(filename = 'F29880.NTR',plotyn = 1, tslice = 0.2);

% Open the file
fid = fopen (filename, 'w');


% Write the header
for k = 1:7
  fdisp(fid, UFILE.header{k});
end

writedata (DATA = UFILE.t, DATASIZE = UFILE.M, number_of_columns = 6, fid);
writedata (DATA = UFILE.data, DATASIZE = UFILE.data_size, number_of_columns = 6, fid);

% Write the bottomer
fdisp(fid, UFILE.bottomer);

% Close the file
fclose(fid); 




function writedata (DATA, DATASIZE, number_of_columns, fid);

  number_of_columns = 6;
  resto = mod(DATASIZE, number_of_columns);
  if (resto == 0)
    fprintf(fid, '  %e %e %e %e %e %e\n', DATA);
  else
    number_of_rows = floor(DATASIZE/number_of_columns) + 1;
    k_start = 1;
    k_end = k_start - 1 + number_of_columns;
    for k = 1:number_of_rows - 1
      fprintf(fid, '  %e %e %e %e %e %e\n', DATA(k_start:k_end));  
      k_start = k*number_of_columns + 1;
      k_end = (k+1)*number_of_columns;
    end
    k_start = k*number_of_columns + 1;
    k_end = k_start - 1 + resto;
    fprintf(fid, ' ');
    for k = k_start:k_end
      fprintf(fid, ' %e', DATA(k));  
    end
    fprintf(fid,  '\n');  
  end

endfunction

endfunction

