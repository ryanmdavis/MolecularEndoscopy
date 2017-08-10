% each dataset must be in it's own folder

function save_file_path = processRaman(dir_in,varargin)

% assign optional arguments
% use 'shift_std_spectra_wavenumber',-10 for Renishaw
% if 'input_recon_param_location' == 1 then we are specifying the location
% of the *reconstruction parameters.mat file
invar = struct('use_dB',0,'redo',0,'shift_std_spectra_wavenumber',0,'show_images',0,'num_background_pc',3,'input_rp_file_location',0);
argin = varargin;
invar = generateArgin(invar,argin);

if invar.input_rp_file_location
    load(dir_in,'out');
    parent_dir=out.rp.parent_dir;
    rp=out.rp;
    im1=out.im1;
    im2=out.im2;
else
    parent_dir=dir_in;
    % make sure parent dir has slash at end
    if ~strcmp(parent_dir(end),mkslash) parent_dir = strcat(parent_dir,mkslash); end
    rp=struct('im1_path','','im2_path','','reg_param',[],'im1',[],'im2',[],'im2_nrow',[],'im2_ncol',[],'fused',[],'tissue_map',[],'A',[],'wavenumber',[],'channel_names',[],'parent_dir',parent_dir);
end



% % if the reconstruction parameters for this dataset already exist, load
% % them, if not then define the structure.
% if ~exist(strcat(parent_dir,'reconstruction parameters.mat'),'file') || invar.redo==-1
%     rp=struct('im1_path','','im2_path','','reg_param',[],'im1',[],'im2',[],'im2_nrow',[],'im2_ncol',[],'fused',[],'tissue_map',[],'A',[],'wavenumber',[],'channel_names',[]);
% else
%     rp = loadRp(parent_dir);
% end

% load Raman data
if isempty(rp.im2_path) || invar.redo==3 || invar.redo==-1
    [filename,pathname] = uigetfile({'*.txt','Endoscope data';'*.*','All Files';
            '*.spc','Microscope data';},'Locate Image 2',parent_dir);
    % make sure parent dir has slash at end
    if ~strcmp(pathname(end),mkslash) pathname = strcat(pathname,mkslash); end
    rp.im2_path=pathname;
    rp.im2_filename=filename;
    dot_loc=strfind(rp.im2_filename,'.');
    filename_save=strcat(rp.im2_path,rp.im2_filename(1:dot_loc(end)-1));
    
    % save recon parameters
    saveRp(filename_save,rp);
else
    % define path to save the processed Raman data (.pr.mat)
    dot_loc=strfind(rp.im2_filename,'.');
    filename_save=strcat(rp.im2_path,rp.im2_filename(1:dot_loc(end)-1));
end

% load first image
if isempty(rp.im1_path) || invar.redo==1 || invar.redo==-1
    [filename,pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'Locate Image 1',...
          parent_dir);
    if pathname ~= 0
        parent_dir=pathname;
        rp.im1_path=strcat(pathname,filename);
        im1=imread(rp.im1_path);
    else
        rp.im1_path='none';
        im1=zeros(256);
    end
    saveRp(filename_save,rp);
elseif strcmp(rp.im1_path,'none')
    im1=zeros(256);
end


% draw user specified mask that surrounds tissue
if ~strcmp(rp.im1_path,'none') && (isempty(rp.tissue_map) || invar.redo==2 || invar.redo==-1)
    f_p=figure;
    imagesc(im1);
    title({'Draw ROI surrounding tissue','but excluding fiducials'});
    axis image off
    xlabel('right click and "create mask"');
    niceFigure(f_p);
    rp.tissue_map=roipoly;
    close(f_p);
end




% recon second image based on if this the data were separated into rows
if isempty(strfind(rp.im2_filename,'_to_')) && ~isempty(strfind(rp.im2_filename,'.spc'))
    spectra=readRenishawSpc(strcat(rp.im2_path,rp.im2_filename));
    rows_cols=[];
    system='Renishaw';
elseif ~isempty(strfind(rp.im2_filename,'.spc'))
    [spectra,rows_cols]=concatenateRamanSpectraRows(rp.im2_path);
    rp.im2_nrow=rows_cols(1);
    rp.im2_ncol=rows_cols(2);
    system='Renishaw';
elseif ~isempty(strfind(rp.im2_filename,'.txt'))
    spectra = readRamanEndoscopeTxt(strcat(rp.im2_path,rp.im2_filename));
    rows_cols=[];
    system='Endoscope';
end

% get info about image matrix size
if isempty(rows_cols) && (isempty(rp.im2_nrow) || invar.redo==4 || invar.redo==-1)
    display(strcat('Total spectra: ',num2str(size(spectra.spectra,1))));
    rp.im2_nrow=input('Enter the number of rows in the Raman map: ');
    saveRp(filename_save,rp);
end

if isempty(rows_cols) && (isempty(rp.im2_ncol) || invar.redo==4 || invar.redo==-1)
    rp.im2_ncol=input('Enter the number of cols in the Raman map: ');
    saveRp(filename_save,rp);
end

% reshape and arrange data for unmixing, unless raw data was split up into
% rows 
if isempty(rows_cols)
    b=permute(reshape(spectra.spectra,rp.im2_ncol,rp.im2_nrow,size(spectra.spectra,2)),[2,1,3]);
else
    b=spectra.spectra;
end

% if first image (i.e. photograph) is not loaded, set tissue mask as
% "everything"
if strcmp(rp.im1_path,'none')
   rp.tissue_map=ones(rp.im2_nrow,rp.im2_ncol);
end

%% define forward problem A (b=Ax) and solve for x
if isempty(rp.A) || invar.redo==6 || invar.redo==-1
    % [A,channel_names,wavenumber]=defineRamanForwardProblem(system,invar.num_background_pc);
    [A,wavenumber,channel_names,num_np_channels]=queryA;
    
    % if specified, shift standard spectra x axis
    wn_per_bin=wavenumber(2)-wavenumber(1);
    for channel_num=1:size(A,1)%-size(background_pc,1)
        A(channel_num,:)=interp1(wavenumber,A(channel_num,:),wavenumber+invar.shift_std_spectra_wavenumber/wn_per_bin);
    end
    A(isnan(A))=0; %interp1 puts NaN into blank spots
    rp.A=A;
    rp.wavenumber=wavenumber;
    rp.channel_names=channel_names;
    rp.num_np_channels=num_np_channels;
    rp.concentrations=ones(1,num_np_channels);
    saveRp(filename_save,rp);
end
    


% % remove the edges of the standard spectra (A) since they sometimes drop to
% % zero and interfere with fitting
% spectra_truncation=15;
% rp.A=rp.A(:,1+spectra_truncation:end-spectra_truncation);
% A_wavenumber=rp.wavenumber(:,1+spectra_truncation:end-spectra_truncation);
% out_unmix=pivUnmixing(b,rp.A,spectra.wavenumber,A_wavenumber);

% kludge fix: forthis dataset, one of the rows was shifted!?
if strcmp(parent_dir,'C:\Users\rdavis5\Documents\Gambhir lab\Gambhir Data and Analysis\Raman Endoscope\02172017 - CA9 and CD47 cells\')
    b(22,:,:)=circshift(b(22,:,:),[0 -9 0]);
end

% perform unmixing
out_unmix=pivUnmixing(b,rp.A,spectra.wavenumber,rp.wavenumber);
% if we are specifying different concentrations for each channel
if invar.redo==5 % if we are specifying different concentrations for each channel
    rp.concentrations=queryRelativeConcGUI(rp.channel_names(1:rp.num_np_channels),rp.concentrations);
end

%% perform registration
% find out the scale,shift, and rotation needed to match the two images
if ~strcmp(rp.im1_path,'none') && (isempty(rp.reg_param) || invar.redo==6 || invar.redo==-1)
    rp.reg_param = determineAffineTransform2(im1,squeeze(out_unmix.x(1,:,:)));
    saveRp(filename_save,rp);
elseif isempty(rp.reg_param)
    % make im1 a RGB image of zeros of same size as im2
    im1 = zeros([size(squeeze(out_unmix.x(1,:,:))) 3]);
    rp.reg_param=struct('shift_rc',[0 0],'rotation',0,'scale',1,'backwards_rot',@(x,y) round([y x]),'center_im1',round([size(im1,1) size(im1,2)]),'padded_im_size_rc',[size(im1,1) size(im1,2)],'padded_im_center_xy',round([size(im1,2) size(im1,1)]/2),'im2_pad_size_rc',[size(im1,1) size(im1,2)]);
    saveRp(filename_save,rp);
end

% for each channel, do the registration and overlay
mf=mfilename('fullpath');
slash_loc=strfind(mf,mkslash);
colormap_path=strcat(mf(1:slash_loc(end)),'colormap',mkslash,'my colormaps.mat');
load(colormap_path,'colormap_bg');  %load the colormap

% out(1:size(out_unmix.x,1))=struct('fused',[],'im1_reg',[],'im2_reg',[],'im1',im1,'im2',out_unmix,'scale',[],'colormap',colormap_bg,'reg_param',rp.reg_param,'A',A,'x',out_unmix.x);
im_reg_blank=zeros([size(out_unmix.x,1) rp.reg_param.padded_im_size_rc 3]);
im_reg_blank_grayscale=zeros([size(out_unmix.x,1) rp.reg_param.padded_im_size_rc]);

% define structure
out=struct('fused',im_reg_blank,'im1_reg_rgb',im_reg_blank,'im2_reg_rgb',im_reg_blank,'im2_reg_grayscale',im_reg_blank_grayscale,'im1',im1,'im2',out_unmix,'scale',zeros(size(out_unmix.x,1),2),'colormap',colormap_bg,'reg_param',rp.reg_param,'A',rp.A,'x',out_unmix.x,'channel_names',{rp.channel_names},'use_dB',invar.use_dB,'pr_file_loc',[],'rp',rp);
for channel_num=1:size(out_unmix.x,1)
    % determine standarding range for Raman map
%     [im1_out, im2_out]=registerImagesByShiftingGrayScale(rp.tissue_map,squeeze(out_unmix.x(1,:,:)),rp.reg_param);
    [reg_mask_out, im2_out]=registerImagesByShifting2(rp.tissue_map,squeeze(out_unmix.x(1,:,:)),rp.reg_param);
    im2_out(im2_out<0)=1e-10;
    masked_im2=im2_out.*reg_mask_out; % apply ROI mask
    if invar.use_dB
        max_raman_intensity=20*log10(max(max(masked_im2)));
        out.scale(channel_num,:)=[max_raman_intensity-60 max_raman_intensity];
        channel_image=squeeze(out_unmix.x(channel_num,:,:));
        channel_image(channel_image<0)=1e-10;
        raman_grayscale=20*log10(channel_image);
        raman_rgb=intensity2RGB(raman_grayscale,colormap_bg,out.scale(channel_num,:));
    else
        max_raman_intensity=max(max(masked_im2));
        out.scale(channel_num,:)=[max_raman_intensity*0.05 max_raman_intensity*0.95];
        raman_grayscale=squeeze(out_unmix.x(channel_num,:,:));
        raman_rgb=intensity2RGB(raman_grayscale,colormap_bg,out.scale(channel_num,:));
    end

    % register and write output.  Show image only once if specified
    if channel_num==1 && invar.show_images
        show_boolean=1;
    else
        show_boolean=0;
    end
    [im1_out_rgb, im2_out_rgb]=registerImagesByShifting2(im1,raman_rgb,rp.reg_param,'show_images',show_boolean);
    out.fused(channel_num,:,:,:) = double(imfuse(im1_out_rgb,im2_out_rgb,'blend'))/255;
    out.im1_reg_rgb(channel_num,:,:,:) = im1_out_rgb;
    out.im2_reg_rgb(channel_num,:,:,:) = im2_out_rgb;
    
    % register and write output
    [~, im2_out_grayscale]=registerImagesByShifting2(im1,raman_grayscale,rp.reg_param);
    out.im2_reg_grayscale(channel_num,:,:) = reshape(im2_out_grayscale,1,size(im2_out_grayscale,1),size(im2_out_grayscale,2));
end

% generate file name for saving and save
slash_loc=strfind(rp.im2_path,mkslash);
dot_loc=strfind(rp.im2_filename,'.');
filename_save=rp.im2_filename(1:dot_loc(end)-1);
save_file_path=strcat(rp.im2_path,filename_save,'.pr.mat');
out.pr_file_loc=save_file_path;
save(out.pr_file_loc,'out','-v7.3');

% display(strcat('processed data written to: ',save_file_path))
end

function saveRp(filename_save,rp)
    save(filename_save,'rp');
end

function rp = loadRp(filename_save)
    load(filename_save,'rp');
end