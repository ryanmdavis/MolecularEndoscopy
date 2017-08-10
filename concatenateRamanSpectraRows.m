function [spectra_out,rows_cols]=concatenateRamanSpectraRows(parent_dir)

d=dir(parent_dir);
d=d(3:end); % remove . and .. directories

% make sure parent dir has slash at end
if ~strcmp(parent_dir(end),mkslash) parent_dir = strcat(parent_dir,mkslash); end

% load first spectra from slides
spectrum1=readRenishawSpc(strcat(parent_dir,d(1).name));

% find number of .spc files
spc_present=~cellfun('isempty',strfind({d.name},'.spc'));
spc_indicies=1:size(d,1);
spc_indicies=spc_indicies(spc_present);

% find the row number of each file
for file_num=1:size(spc_indicies,2)    
   to_loc=strfind(d(file_num).name,'_to_'); 
   underscore_loc=strfind(d(file_num).name,'_');
   column_number(file_num)=str2double(d(file_num).name(underscore_loc(1)+1:to_loc-1));
end
file_num_to_row_num=column_number/size(spectrum1.spectra,1)+1;

% allocate memory for large dataset
raman_map=zeros(size(spc_indicies,2), size(spectrum1.spectra,1), size(spectrum1.spectra,2));

for file_num=1:size(spc_indicies,2)    
   % load spectrum and unmix
    spectra=readRenishawSpc(strcat(parent_dir,d(file_num).name));
    raman_map(file_num_to_row_num(file_num),:,:)=spectra.spectra;
end

spectra_out=spectrum1;
spectra_out.spectra=raman_map;
rows_cols=[size(spc_indicies,2) size(spectrum1.spectra,1)];