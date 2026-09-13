% Image output test

image_dir = strcat(dir,'images\iamges_',day_for_filename,'_',time_for_filename);
mkdir(image_dir);
cellspos=or(CKpos,CD45pos);
if tub_yesno~=1
    cellspos=or(cellspos,TUBpos);
end
if ar_yesno~=1
    cellspos=or(cellspos,ARpos);
end
for position_step=position_start:position_skip:position_end,

    for zstep=3:slice_skip:z_stack_depth

        cells_in_this_spot = and(locationarray(:,1)==position_step,round(locationarray(:,4))==zstep);
        if dapi_yesno~=1
            cells_in_image=cells_in_this_spot;
        else
            cells_in_image=and(cells_in_this_spot,cellspos);
        end
        if any(cells_in_image)
            % Read in DAPI
            if location==1
                generic_marker = DAPImarker;
                generic_channel = DAPIchannel;
            elseif location==2
                multichannel_generic = multichannelDAPI;
                generic_channelTEXT = DAPIchannelTEXT;
            end
            CTCQ_read_in_for_image;
            DAPI_in_image = double(generic_in);
                       
            % Read in CK
            if location==1
                generic_marker = CKmarker;
                generic_channel = CKchannel;
            elseif location==2
                multichannel_generic = multichannelCK;
                generic_channelTEXT = CKchannelTEXT;
            end
            CTCQ_read_in_for_image;
            CK_in_image = double(generic_in);
            
            % Read in CD45
            if location==1
                generic_marker = CD45marker;
                generic_channel = CD45channel;
            elseif location==2
                multichannel_generic = multichannelCD45;
                generic_channelTEXT = CD45channelTEXT;
            end
            CTCQ_read_in_for_image;
            CD45_in_image = double(generic_in);
            
            % Read in TUB
            %if tub_yesno~=1
            %if location==1
             %   generic_marker = TUBmarker;
              %  generic_channel = TUBchannel;
            %elseif location==2
            %    multichannel_generic = multichannelTUB;
             %   generic_channelTEXT = TUBchannelTEXT;
            %end
            %CTCQ_read_in_for_image;
            %TUB_in_image = double(generic_in);
            %end          
           
            % Read in AR
     %       if ar_yesno~=1
      %      if location==1
       %         generic_marker = ARmarker;
        %        generic_channel = ARchannel;
         %   elseif location==2
          %      multichannel_generic = multichannelAR;
           %     generic_channelTEXT = ARchannelTEXT;
            %end
            %CTCQ_read_in_for_image;
            %AR_in_image = double(generic_in);
            %end      
            
            
            figure(995)
            imshow(DAPI_in_image,[]);
            hold on
           %  DAPI image
                  if dapi_yesno~=1
            plot(locationarray(cells_in_this_spot,2),locationarray(cells_in_this_spot,3),'gs','MarkerSize',8)
            end 
            plot(locationarray(and(cells_in_this_spot,CKpos),2) ,locationarray(and(cells_in_this_spot,CKpos),3) ,'rs','MarkerSize',10)
            plot(locationarray(and(cells_in_this_spot,CD45pos),2) ,locationarray(and(cells_in_this_spot,CD45pos),3) ,'cs','MarkerSize',12)
             if tub_yesno~=1
            plot(locationarray(and(cells_in_this_spot,TUBpos),2) ,locationarray(and(cells_in_this_spot,TUBpos),3) ,'ms','MarkerSize',14)
            end
            if ar_yesno~=1
            plot(locationarray(and(cells_in_this_spot,ARpos),2) ,locationarray(and(cells_in_this_spot,ARpos),3) ,'ys','MarkerSize',16)
            end
             image_filename = strcat(image_dir,'\Labeled_',DAPImarker,'_image_',day_for_filename,'_',time_for_filename,'_z=',num2str(zstep),'_Loc=',num2str(position_step));
            print('-djpeg','-f995',image_filename);
            close 995
            
             figure(995)
            
            imshow(CK_in_image,[]);
            hold on
            % CK image
            if dapi_yesno~=1
            plot(locationarray(cells_in_this_spot,2) ,locationarray(cells_in_this_spot,3) ,'gs','MarkerSize',8)
            end 
            plot(locationarray(and(cells_in_this_spot,CKpos),2) ,locationarray(and(cells_in_this_spot,CKpos),3) ,'rs','MarkerSize',10)
            plot(locationarray(and(cells_in_this_spot,CD45pos),2) ,locationarray(and(cells_in_this_spot,CD45pos),3) ,'cs','MarkerSize',12)
             if tub_yesno~=1
            plot(locationarray(and(cells_in_this_spot,TUBpos),2) ,locationarray(and(cells_in_this_spot,TUBpos),3) ,'ms','MarkerSize',14)
            end
            if ar_yesno~=1
            plot(locationarray(and(cells_in_this_spot,ARpos),2) ,locationarray(and(cells_in_this_spot,ARpos),3) ,'ys','MarkerSize',16)
            end
            
            image_filename = strcat(image_dir,'\Labeled_',CKmarker,'_image_',day_for_filename,'_',time_for_filename,'_z=',num2str(zstep),'_Loc=',num2str(position_step));
            print('-djpeg','-f995',image_filename);
            close 995
            
                     figure(995)
            imshow(CD45_in_image,[]);
            hold on
            % CD45 image
            if dapi_yesno~=1
            plot(locationarray(cells_in_this_spot,2) ,locationarray(cells_in_this_spot,3) ,'gs','MarkerSize',8)
            end
            plot(locationarray(and(cells_in_this_spot,CKpos),2) ,locationarray(and(cells_in_this_spot,CKpos),3) ,'rs','MarkerSize',10)
            plot(locationarray(and(cells_in_this_spot,CD45pos),2) ,locationarray(and(cells_in_this_spot,CD45pos),3) ,'cs','MarkerSize',12)
              if tub_yesno~=1
           plot(locationarray(and(cells_in_this_spot,TUBpos),2) ,locationarray(and(cells_in_this_spot,TUBpos),3) , 'ms','MarkerSize',14)
            end
            if ar_yesno~=1
            plot(locationarray(and(cells_in_this_spot,ARpos),2) ,locationarray(and(cells_in_this_spot,ARpos),3) ,'ys','MarkerSize',16)
            end
            image_filename = strcat(image_dir,'\Labeled_',CD45marker,'_image_',day_for_filename,'_',time_for_filename,'_z=',num2str(zstep),'_Loc=',num2str(position_step));
            print('-djpeg','-f995',image_filename);
            close 995
            
      %       if tub_yesno~=1
      %     figure(995)
      %     imshow(TUB_in_image,[]);
      %      hold on
            %  TUB image
      %      plot(locationarray(cells_in_this_spot,2) ,locationarray(cells_in_this_spot,3) ,'g*')
      %      plot(locationarray(and(cells_in_this_spot,CKpos),2) ,locationarray(and(cells_in_this_spot,CKpos),3) ,'rs')
       %     plot(locationarray(and(cells_in_this_spot,CD45pos),2) ,locationarray(and(cells_in_this_spot,CD45pos),3) ,'r+')
       %     plot(locationarray(and(cells_in_this_spot,TUBpos),2) ,locationarray(and(cells_in_this_spot,TUBpos),3) ,'m.')
       %     if ar_yesno~=1
        %    plot(locationarray(and(cells_in_this_spot,ARpos),2) ,locationarray(and(cells_in_this_spot,ARpos),3) ,'yd')
        %    end
        %    image_filename = strcat(image_dir,'\Labeled_',TUBmarker,'_image_',day_for_filename,'_',time_for_filename,'_z=',num2str(zstep),'_Loc=',num2str(position_step));
        %    print('-djpeg','-f995',image_filename);
        %    close 995
        %     end
            
       %       if ar_yesno~=1
       %              figure(995)
       %     imshow(AR_in_image,[]);
       %     hold on
             % AR image
       %     plot(locationarray(cells_in_this_spot,2) ,locationarray(cells_in_this_spot,3) ,'g*')
       %     plot(locationarray(and(cells_in_this_spot,CKpos),2) ,locationarray(and(cells_in_this_spot,CKpos),3) ,'rs')
       %     plot(locationarray(and(cells_in_this_spot,CD45pos),2) ,locationarray(and(cells_in_this_spot,CD45pos),3) ,'r+')
       %      if tub_yesno~=1
       %     plot(locationarray(and(cells_in_this_spot,TUBpos),2) ,locationarray(and(cells_in_this_spot,TUBpos),3) ,'m.')
       %      end
       %     plot(locationarray(and(cells_in_this_spot,ARpos),2) ,locationarray(and(cells_in_this_spot,ARpos),3) ,'yd')

       %     image_filename = strcat(image_dir,'\Labeled_',ARmarker,'_image_',day_for_filename,'_',time_for_filename,'_z=',num2str(zstep),'_Loc=',num2str(position_step));
       %     print('-djpeg','-f995',image_filename);
          
 
       %     close 995
        %      end
        end
        
    end     
   
end;