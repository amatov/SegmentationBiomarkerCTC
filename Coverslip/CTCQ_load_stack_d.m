
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DAPI Z STACK AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND DAPI THRESHOLD OF Z STACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOSAIC IN GEDI CHIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _c: generalized to load any type of stack as long as the parameters are
% defined before this script is called
% moved all of the actual loading and filtering to CTCQ_LogFilterThresh 
% this script just opens the files
% _b: changed "level" terminology to level_DAPI to reflect thresholding for
% all channels in cell definition
% _a: pre-filters with difference of Gaussians

slice_counter=0;% this tells the program where to put each slice.  used when skipping slices.


generic_zstack = [];
generic_zstackbw = [];

if location==1,% if WCMC
    for zstep=1:slice_skip:z_stack_depth,% loop through z stack, load images, make stack
        slice_counter=slice_counter+1;
        fprintf(1,strcat(['Loading ' generic_marker ' slice %i of %i of position %i\n']),zstep,z_stack_depth,position_step);
        generic_name=strcat(dapistart,'p000',pindex(position_step,:),'t00000001z0',zindex(zstep,:),generic_channel,'.TIF');
        generic_file=strcat(dir,'\',generic_name);
        generic_in=imread(generic_file,'tiff');%read dapi image
        
       

%%subtract background
generic=double(generic_in);

%%put images in stack
generic_zstack(:,:,slice_counter)=generic;
       [detectionResults,qq_thresholded] =spotDetector(generic);
       generic_zstackbw(:,:,slice_counter)=qq_thresholded;
        %CTCQ_LogFilterThresh;
    end

elseif location==2,% if Cornell Ithaca
   generic_file=strcat([dir '\' dapistart '_' generic_channelTEXT '_L' num2str(position_step) '_Sum.lsm']);
   generic_temp=tiffread(generic_file);
   % if there are multiple channels in this file, set generic_temp

   
   
    for zstep=1:slice_skip:z_stack_depth,% loop through z stack, load images, make stack
        slice_counter = slice_counter+1; 

        if multichannel_generic==0 %generic_temp is different depending on whether the input image has 1 or multiple channels
            generic_in = generic_temp(zstep).data;%pull out DAPI image from lsm stack 
        else
            generic_in = generic_temp(zstep).data{multichannel_generic};
        end
        
        log_generic_in = log(double(generic_in)+1);

%%subtract background
generic=double(generic_in);

%%put images in stack
generic_zstack(:,:,slice_counter)=generic;
       [detectionResults,qq_thresholded] =spotDetector(generic);
       generic_zstackbw(:,:,slice_counter)=qq_thresholded;
       % CTCQ_LogFilterThresh;
    end
end;


