function [im2,rp]=loadRawRamanFromHD(rp)
    % recon second image if this the data were separated into rows
    if isempty(strfind(rp.im2_filename,'_to_')) && ~isempty(strfind(rp.im2_filename,'.spc'))
        im2=readRenishawSpc(strcat(rp.im2_path,rp.im2_filename));
        im2 = queryAndReshape(im2,rp);
    elseif ~isempty(strfind(rp.im2_filename,'.spc'))
        f_h=figure;
        imagesc(ones(1,10));
        axis image
        text(0,0,'Loading raw Raman data from hard drive.  Please wait.');
        [im2,rows_cols]=concatenateRamanSpectraRows(rp.im2_path);
        rp.im2_nrow=rows_cols(1);
        rp.im2_ncol=rows_cols(2);
        close(f_h)
    elseif ~isempty(strfind(rp.im2_filename,'.txt'))
        im2 = readRamanEndoscopeTxt(strcat(rp.im2_path,rp.im2_filename));
        im2 = queryAndReshape(im2,rp);
        rp.im2_nrow=size(im2.spectra,1);
        rp.im2_ncol=size(im2.spectra,2);        
    else
        rp.im2_nrow=size(im2.b,1); 
        rp.im2_ncol=size(im2.b,2);
    end
end


function im2 = queryAndReshape(im2,rp)
    % get info about image matrix size
    if isempty(rp.tissue_map)
        if isempty(rp.im2_nrow) %|| invar.redo==4 || invar.redo==-1
            display(strcat('Total spectra: ',num2str(size(im2.spectra,1))));
            rp.im2_nrow=input('Enter the number of rows in the Raman map: ');
        end
        if isempty(rp.im2_ncol) %|| invar.redo==4 || invar.redo==-1
            rp.im2_ncol=input('Enter the number of cols in the Raman map: ');
        end
        im2.spectra = permute(reshape(im2.spectra,rp.im2_ncol,rp.im2_nrow,size(im2.spectra,2)),[2,1,3]);
    else %if we've already loaded this dataset before, the image size is stored under rp.tissue_map
        im2.spectra = permute(reshape(im2.spectra,size(rp.tissue_map,1),size(rp.tissue_map,2),size(im2.spectra,2)),[2,1,3]);
    end
end