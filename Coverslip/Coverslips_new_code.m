clear all
%[filename2,dir2] = uigetfile('*.xls' , 'Chose the excel file that you want to re-generate');
%filename=strcat(dir2,filename2);
%coordinates=xlsread(filename);
coordinates1=xlsread('C:\Users\HRISTOV\Desktop\Pt.9283-920_Cover1\Guang 1-90 ctc positions - coverslips.xlsx');
all_coordinates=[];
all_info=[];

for i=24:24
    meanCTCintensity=[];
    ctc_location=[];
    meanCTCarea=[];
    coordinates=[];
    ctc1_location=[];
    info=[];
DAPI_in_image= imread(strcat('C:\Users\HRISTOV\Desktop\Pt.9283-920_Cover1\Pt.9283-920_Cover1_z4m',sprintf('%03d',i),'c4.tif'),'tiff');

CTC_in_image= imread(strcat('C:\Users\HRISTOV\Desktop\Pt.9283-920_Cover1\Pt.9283-920_Cover1_z4m',sprintf('%03d',i),'c2.tif'),'tiff');

%TUB_in_image= imread('C:\Users\HRISTOV\Desktop\Pt.9283-920_Cover1\Pt.9283-920_Cover1_z4m096c3.tif','tiff');
 
%AR_in_image= imread('C:\Users\HRISTOV\Desktop\ChipA0303-0001.tif_Files\ChipA0303-0001_p000002t00000001z002c03.tif','tiff');

L_in_image= imread(strcat('C:\Users\HRISTOV\Desktop\Pt.9283-920_Cover1\Pt.9283-920_Cover1_z4m',sprintf('%03d',i),'c5.tif'),'tiff');

wavelet_image=double(DAPI_in_image);

wavelet_CTC_image=double(CTC_in_image);

wavelet_L_image=double(L_in_image);

%wavelet_TUB_image=double(TUB_in_image);

%wavelet_AR_image=double(AR_in_image);
%cutI=wavelet_image;
% [cutoffInd, cutDAPI] = cutFirstHistMode(cutI,0);
% qq_thresholded = cutI>cutDAPI*2;

[detectionResults,qq_thresholded] =spotDetector(wavelet_image);
 DAPIzstackbw(:,:,:)=qq_thresholded;
 mask = bwlabel(DAPIzstackbw);
 DAPIzstackstats = regionprops(mask, DAPIzstackbw, {'Area','Centroid','MeanIntensity'});
[num_cells junkthree]=size(DAPIzstackstats);

         


for junk=1:num_cells,
    

DAPIarea(junk)=DAPIzstackstats(junk).Area;
DAPIintensity(junk)=DAPIzstackstats(junk).MeanIntensity;

locationarray(junk,:) = [2 DAPIzstackstats(junk).Centroid(:)'];
end

cutI=DAPIintensity.*DAPIarea;
[cutoffInd, cutDAPI] = cutFirstHistMode(cutI,0);


a=[DAPIzstackstats.MeanIntensity].*[DAPIzstackstats.Area];

 dapi_selected = ismember(mask, find([a]>cutDAPI*0.57 & [DAPIzstackstats.Area] < 40));
 mask=bwlabel(dapi_selected);


%[detectionResults,tub_thresholded] =spotDetector(wavelet_TUB_image);
% TUBzstackbw(:,:,:)=tub_thresholded;
% mask_tub = logical(TUBzstackbw);
% TUBzstackstats = regionprops(mask_tub, TUBzstackbw, {'Area','Centroid','MeanIntensity'});
%[num_tub junkthree]=size(TUBzstackstats);


%[detectionResults,ar_thresholded] =spotDetector(wavelet_AR_image);
% ARzstackbw(:,:,:)=ar_thresholded;
% mask_ar = logical(ARzstackbw);
% ARzstackstats = regionprops(mask_ar, ARzstackbw, {'Area','Centroid','MeanIntensity'});
%[num_ar junkthree]=size(ARzstackstats);

%cutI=wavelet_L_image;
%[cutoffInd, cutDAPI] = cutFirstHistMode(cutI,0);
%l_thresholded = cutI>cutDAPI*1.2;


 [detectionResults,l_thresholded] =spotDetector(wavelet_L_image);
 Lzstackbw(:,:,:)=l_thresholded;
 mask_l = bwlabel(Lzstackbw);
 Lzstackstats = regionprops(mask, l_thresholded, {'MeanIntensity','Area','FilledArea'});

 
 [detectionResults,ctc_thresholded] =spotDetector(wavelet_CTC_image);
 %[cutoffInd, cutCTC] = cutFirstHistMode(wavelet_CTC_image,0);
 % ctc_thresholded = wavelet_CTC_image>cutCTC*1.5;
 CTCzstackbw(:,:,:)=ctc_thresholded;
 mask_ctc = bwlabel(CTCzstackbw);
 
   %indices = find(ctc_thresholded()==0);
    %      ctc_thresholded(indices) = [];
%coeff=var(ctc_thresholded(:));
 
 %mask_new=mask_ctc;
%mask_new = and(mask,mask_ctc);
% mask_new =  and(mask_new,mask_tub);
 % mask_new =  and(mask_new,mask_ar);
 % mask_l=and(mask,mask_l);
  %new_mask=xor(mask_new,mask_l);
  %new_mask=and(new_mask,mask_new);
  
 % cutI=CTCzstackbw*mask_new;
 % [cutoffInd, cutCTC] = cutFirstHistMode(cutI,0);
 %CTCzstackbw(:,:,:) = cutI>cutCTC*0;
 % mask_new=logical(CTCzstackbw);
  
  
 CTCzstackstats = regionprops(mask, CTCzstackbw, {'MeanIntensity','Area','Centroid','FilledArea'});
% DAPIstats = regionprops(mask_new, DAPIzstackbw, {'Area','Centroid','MeanIntensity'});
% TUBstats = regionprops(mask_new, TUBzstackbw, {'Area','Centroid','MeanIntensity'});
% ARstats = regionprops(mask_new, ARzstackbw, {'Area','Centroid','MeanIntensity'});
% CTCstats = regionprops(mask_new, CTCzstackbw, {'Area','Centroid','MeanIntensity'});
 [num_ctc junkthree]=size(CTCzstackstats);

 

 






cell_counter=0;
for junk=1:num_cells,
DAPIarea(junk)=DAPIzstackstats(junk).Area;
DAPIintensity(junk)=DAPIzstackstats(junk).MeanIntensity;
locationarray(junk,:) = [2 DAPIzstackstats(junk).Centroid(:)'];
end
%for junk=1:num_cells,
%   if area(junk)>16.5 && area(junk)<10000000
%    cell_counter=cell_counter+1;
%    locationarray(cell_counter,:) = [2 DAPIzstackstats(junk).Centroid(:)'];

%            meanLarray(cell_counter) = Lzstackstats(junk).MeanIntensity;
%          meanLarea(cell_counter) = Lzstackstats(junk).Area;
% end
%end
cell_counter=0;
 
for junk=1:num_ctc
    
    cell_counter=cell_counter+1;
   % ctc_location(cell_counter,:) = [i CTCzstackstats(junk).Centroid(:)'];
    %       meanCTCintensity(cell_counter) = CTCzstackstats(junk).MeanIntensity;
     %      meanCTCarea(cell_counter) = CTCzstackstats(junk).Area;
           
      meanCD45array(cell_counter) =  Lzstackstats(junk).MeanIntensity;
       meanCTCarray(cell_counter) = CTCzstackstats(junk).MeanIntensity;    
           locationarray1(cell_counter,:) = [2 DAPIzstackstats(junk).Centroid(:)'];
           
            CD45area(cell_counter) =  Lzstackstats(junk).Area;
            CTCarea(cell_counter) = CTCzstackstats(junk).Area;    
     
             CD45fill(cell_counter) =  Lzstackstats(junk).FilledArea;
            CTCareafill(cell_counter) = CTCzstackstats(junk).FilledArea;   
           %DAPIintensity(cell_counter) = DAPIzstackstats(junk).MeanIntensity;
           %DAPIarea(cell_counter) = DAPIzstackstats(junk).Area;
           
           %TUBintensity(cell_counter) = TUBstats(junk).MeanIntensity;
           %TUBarea(cell_counter) = TUBstats(junk).Area;
           
           %ARintensity(cell_counter) = ARstats(junk).MeanIntensity;
           %ARarea(cell_counter) = ARstats(junk).Area;
           
          % CTCintensity(cell_counter) = DAPIstats(junk).MeanIntensity;
          % CTCarea(cell_counter) = DAPIstats(junk).Area;
           
    
end



%info=[ctc_location meanCTCintensity'];
%info=[info meanCTCarea'];
%all_info=[all_info;info];
% cutI=DAPIintensity.*DAPIarea;
% [cutoffInd, cutDAPI] = cutFirstHistMode(cutI,0);
% DapiCUT = cutI>cutDAPI*0.1;

%if (num_ctc>0)    
 %   if (num_ctc==1)
  %      CtcCUT=1;
   %     meanCTCintensityCUT=1;
   % else

cutI1=meanCTCarray;
[cutoffInd, cutCTC] = cutFirstHistMode(cutI1,0);
CtcCUT = cutI1>cutCTC*0.8%*(var(cutI1)/old_cutI1);

%old_cutI1=var(cutI1);
cutI=meanCD45array;
[cutoffInd, cutCTC] = cutFirstHistMode(cutI,0);%
CD45CUT = cutI>cutCTC*0.5;
br=0;
for i=1:num_ctc;
    if and(CtcCUT(i)==1,CD45CUT(i)==0)
        br=br+1;
        Ctcposition(br)=i;
    end
end
%cutI=meanCTCintensity;
%[cutoffInd, cutCTC] = cutFirstHistMode(cutI,0);
%meanCTCintensityCUT = cutI>cutCTC*100%*(var(cutI)/old_cutI);


%old_cutI=var(cutI);
%end
%CUT=and(meanCTCareaCUT,meanCTCintensityCUT);
    
%totalintensity=(meanCTCintensity+CTCintensity+DAPIintensity+TUBintensity+ARintensity);


 % cutI=totalintensity;
 %[cutoffInd, cutCTC] = cutFirstHistMode(cutI,0);
%totalintensityCUT = cutI>cutCTC*0.5;

%totalCUT=and(CUT,totalintensityCUT);
 
 %figure(995)
           % imshow( DAPI_in_image,[]);
            %imshow( DAPI_in_image,[]);
          %  hold on
         %  plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
        %  for i=1:size(locationarray())
             % if and(locationarray(i,3)>2.9,locationarray(i,3)<219)
       %    plot(locationarray(i,2),locationarray(i,3),'gs','MarkerSize',8);
             % end
      %    end
           
   
     %  cd45=0;
     
      ctc1=0;

           for j=1:num_ctc
            if and(meanCTCintensityCUT(j)>0,CtcCUT(j)>0)
                    ctc1=ctc1+1;
                     ctc1_location(ctc1,1)=i;
                   ctc1_location(ctc1,2)=ctc_location(j,2);
                   ctc1_location(ctc1,3)=ctc_location(j,3);
              
                 
       end
           end
           if (ctc1>0)
figure(i)
           imshow( CTC_in_image,[]);
            %imshow( DAPI_in_image,[]);
            hold on
            
           % indices = find(coordinates(:,3)==0);
           % coordinates(indices,:) = [];
            
           %plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
         for j=1:size(coordinates1())
              if and(coordinates1(j,1)==i,coordinates1(j,4)==1)
                  
           plot(coordinates1(j,2),coordinates1(j,3),'r*','MarkerSize',8)
  
              end
                 end

         %  ctc1=0;
           
          for k=1:ctc1
              %for j=1:num_cells
           %   if totalCUT(i)>0
        %     if and(CTCzstackstats(i).Centroid(1)>(DAPIzstackstats(j).Centroid(1)-0.5),CTCzstackstats(i).Centroid(1)<(DAPIzstackstats(j).Centroid(1))+0.5)
         %      if and(CTCzstackstats(i).Centroid(2)>(DAPIzstackstats(j).Centroid(2)-0.5),CTCzstackstats(i).Centroid(2)<(DAPIzstackstats(j).Centroid(2))+0.5)
               %   if meanCTCarea(i)>mean(meanCTCarea)-1
               %       if meanCTCintensity(i)>650
                  %  if totalintensity(i)>median(totalintensity)
                    
           
                  %  if CtcCUT(i)>0
                  %  if (meanCTCintensityCUT(i)>0)
                 %   if or(and(meanCTCarea(i)==1,meanCTCintensity(i)>2823),or(meanCTCarea(i)==2,or(meanCTCarea(i)==7,or(meanCTCarea(i)==5,meanCTCarea(i)>12))))%6,meanCTCarea(i)> 14)%,)
                    %   ctc1=ctc1+1;
                    plot(ctc1_location(k,2),ctc1_location(k,3),'gs','MarkerSize',8)
                    coordinates(k,:)=ctc1_location(k,:);
                    
             %  text(ctc_location(i,2)+5,ctc_location(i,3)+5,num2str(meanCTCarea(i)),'Color',[0 1 0])
              %     text(ctc_location(i,2)+10,ctc_location(i,3)+10,num2str(meanCTCintensity(i)),'Color',[0 1 1])
                   
                    
               %  end
               % end
              end;
%               figure(1000+i)
%            imshow( L_in_image,[]);
%             %imshow( DAPI_in_image,[]);
%             hold on
%             
%            % indices = find(coordinates(:,3)==0);
%            % coordinates(indices,:) = [];
%             
%            %plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
%    
% 
%          %  ctc1=0;
%            
%           for k=1:ctc1
%               %for j=1:num_cells
%            %   if totalCUT(i)>0
%         %     if and(CTCzstackstats(i).Centroid(1)>(DAPIzstackstats(j).Centroid(1)-0.5),CTCzstackstats(i).Centroid(1)<(DAPIzstackstats(j).Centroid(1))+0.5)
%          %      if and(CTCzstackstats(i).Centroid(2)>(DAPIzstackstats(j).Centroid(2)-0.5),CTCzstackstats(i).Centroid(2)<(DAPIzstackstats(j).Centroid(2))+0.5)
%                %   if meanCTCarea(i)>mean(meanCTCarea)-1
%                %       if meanCTCintensity(i)>650
%                   %  if totalintensity(i)>median(totalintensity)
%                     
%            
%                   %  if CtcCUT(i)>0
%                   %  if (meanCTCintensityCUT(i)>0)
%                  %   if or(and(meanCTCarea(i)==1,meanCTCintensity(i)>2823),or(meanCTCarea(i)==2,or(meanCTCarea(i)==7,or(meanCTCarea(i)==5,meanCTCarea(i)>12))))%6,meanCTCarea(i)> 14)%,)
%                     %   ctc1=ctc1+1;
%                     plot(ctc1_location(k,2),ctc1_location(k,3),'gs','MarkerSize',8)
%                     coordinates(k,:)=ctc1_location(k,:);
%                     
%              %  text(ctc_location(i,2)+5,ctc_location(i,3)+5,num2str(meanCTCarea(i)),'Color',[0 1 0])
%               %     text(ctc_location(i,2)+10,ctc_location(i,3)+10,num2str(meanCTCintensity(i)),'Color',[0 1 1])
%                    
%                     
%                %  end
%                % end
%               end;
%               
%               figure(2000+i)
%            imshow( DAPI_in_image,[]);
%             %imshow( DAPI_in_image,[]);
%             hold on
%             
%            % indices = find(coordinates(:,3)==0);
%            % coordinates(indices,:) = [];
%             
%            %plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
%    
% 
%          %  ctc1=0;
%            
%           for k=1:ctc1
%               %for j=1:num_cells
%            %   if totalCUT(i)>0
%         %     if and(CTCzstackstats(i).Centroid(1)>(DAPIzstackstats(j).Centroid(1)-0.5),CTCzstackstats(i).Centroid(1)<(DAPIzstackstats(j).Centroid(1))+0.5)
%          %      if and(CTCzstackstats(i).Centroid(2)>(DAPIzstackstats(j).Centroid(2)-0.5),CTCzstackstats(i).Centroid(2)<(DAPIzstackstats(j).Centroid(2))+0.5)
%                %   if meanCTCarea(i)>mean(meanCTCarea)-1
%                %       if meanCTCintensity(i)>650
%                   %  if totalintensity(i)>median(totalintensity)
%                     
%            
%                   %  if CtcCUT(i)>0
%                   %  if (meanCTCintensityCUT(i)>0)
%                  %   if or(and(meanCTCarea(i)==1,meanCTCintensity(i)>2823),or(meanCTCarea(i)==2,or(meanCTCarea(i)==7,or(meanCTCarea(i)==5,meanCTCarea(i)>12))))%6,meanCTCarea(i)> 14)%,)
%                     %   ctc1=ctc1+1;
%                     plot(ctc1_location(k,2),ctc1_location(k,3),'gs','MarkerSize',8)
%                     coordinates(k,:)=ctc1_location(k,:);
%                     
%              %  text(ctc_location(i,2)+5,ctc_location(i,3)+5,num2str(meanCTCarea(i)),'Color',[0 1 0])
%               %     text(ctc_location(i,2)+10,ctc_location(i,3)+10,num2str(meanCTCintensity(i)),'Color',[0 1 1])
%                    
%                     
%                %  end
%                % end
%               end;
              
              all_coordinates=[all_coordinates;coordinates];
              
              end;
              
          end;
          
        
          
%end
%              
% 
%           
%        %   ctc=0;
%          figure(997)
%          % imshow( CTC_in_image,[]);
%           
%             imshow( DAPI_in_image,[]);
%             hold on
%          
%           % plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
%             
%           hold on
%           %  plot (and(qq_thresholded,ctc_thresholded),'gs','MarkerSize',8)
%           for i=1:num_cells
%          %     if(meanCTCarray(i)> (max(meanCTCarray)+min(meanCTCarray(meanCTCarray>0)))/2)
%           %        ctc=ctc+1;
%               
% %          if DapiCUT(i)>0
%           plot(locationarray(i,2),locationarray(i,3),'gs','MarkerSize',8)
%           hold on
%           plot (240,193,'rs','MarkerSize',8)
% %          end;
%               
%               
%           end;
%           
%              figure(998)
%             %   imshow( ctc_thresholded,[0 500]);
%             imshow( DAPI_in_image,[]);
%             hold on
%             
%            % indices = find(coordinates(:,3)==0);
%           %  coordinates(indices,:) = [];
%             
%            plot(coordinates(:,1),coordinates(:,2),'r*','MarkerSize',8)
%             
%          %  ctc1=0;
%            
%          % for i=1:num_ctc
%          for i=1:num_cells
%              
%          %     plot(ctc_location(i,2),ctc_location(i,3),'gs','MarkerSize',8)
%                plot(locationarray(i,2),locationarray(i,3),'gs','MarkerSize',8)  
%          end;
%           %     for i=1:size(locationarray())
%              % if and(locationarray(i,3)>2.9,locationarray(i,3)<219)
%           % plot(locationarray(i,2),locationarray(i,3),'ys','MarkerSize',10);
%              % end
%          % end
%         
%         
% 
%             