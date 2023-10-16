% number of sheets
num_sheets = 2^n_LABEL;
[col,row] = meshgrid(1:n_LABEL,1:num_sheets);
pos_neg_matrix = logical(floor(mod(row-1,2.^col)./(2.^(col-1))));

% contains a 1 in the column if the cell should be + for corresponding
% column in pos_neg_data
n_cells = length(meanCD45array);

list_of_namevecs = cell(num_sheets,1);
list_of_cell_counts = zeros(num_sheets,1);        

if (n_LABEL>2)
     pos_neg_matrix=fliplr(pos_neg_matrix);
     export_data(:,length(export_data(1,:)))=[];
end;

for ii=1:num_sheets
    % for each sheet, see which cells match the correct combo or + or - for
    % the corresponding combo of + or - fluorophore
    writecellarray = false(size(meanCD45array));
    for jj=1:n_cells
        writecellarray(jj) = all(pos_matrix(jj,:) == pos_neg_matrix(ii,:));
    end
    cell_subset_locations = export_data(writecellarray,:);
    plusorminus_text = cell(1,n_LABEL);
    plusorminus_text(pos_neg_matrix(ii,:)) = {'+.'}; % excel doesn't like '/' in sheet names
    plusorminus_text(~pos_neg_matrix(ii,:)) = {'-.'};
    together = [label_text' plusorminus_text']';
    reshaped = reshape(together,1,length(label_text)+length(plusorminus_text));
    namevec = strcat(reshaped{:});
    %remove the last / from the sheet title
    namevec = namevec(1:(end-1));
    cell_count_ofthistype = sum(writecellarray);
    cell_count = {strcat(namevec, ' count:') cell_count_ofthistype};
    
       if and(tub_yesno==1,ar_yesno==1)
        header = {'Location Number' 'Z Slice' 'Local Pixel Row' 'Local Pixel Column' 'Global Parallel to Flow (um)' 'Global Perpindicular to Flow (um)' 'Diameter' [CD45marker ' value'] [CKmarker ' value'] 'flag'};
        counts=[counts; cell_count_ofthistype mean(cell_subset_locations(:,7)) std(cell_subset_locations(:,7))  median(cell_subset_locations(:,7)) mean(cell_subset_locations(:,8)) std(cell_subset_locations(:,8))  median(cell_subset_locations(:,8)) mean(cell_subset_locations(:,9)) std(cell_subset_locations(:,9)) median(cell_subset_locations(:,9))];
        elseif (ar_yesno==1)
        header = {'Location Number' 'Z Slice' 'Local Pixel Row' 'Local Pixel Column' 'Global Parallel to Flow (um)' 'Global Perpindicular to Flow (um)' 'Diameter' [CD45marker ' value'] [CKmarker ' value'] [TUBmarker ' Value'] 'flag'};
        counts=[counts; cell_count_ofthistype mean(cell_subset_locations(:,7)) std(cell_subset_locations(:,7))  median(cell_subset_locations(:,7)) mean(cell_subset_locations(:,8)) std(cell_subset_locations(:,8)) median(cell_subset_locations(:,8)) mean(cell_subset_locations(:,9)) std(cell_subset_locations(:,9)) median(cell_subset_locations(:,9)) mean(cell_subset_locations(:,10)) std(cell_subset_locations(:,10)) median(cell_subset_locations(:,10))];
        elseif (tub_yesno==1)
        header = {'Location Number' 'Z Slice' 'Local Pixel Row' 'Local Pixel Column' 'Global Parallel to Flow (um)' 'Global Perpindicular to Flow (um)' 'Diameter' [CD45marker ' value'] [CKmarker ' value'] [ARmarker ' Value'] 'flag'};
        counts=[counts; cell_count_ofthistype mean(cell_subset_locations(:,7)) std(cell_subset_locations(:,7))  median(cell_subset_locations(:,7)) mean(cell_subset_locations(:,8)) std(cell_subset_locations(:,8)) median(cell_subset_locations(:,8)) mean(cell_subset_locations(:,9)) std(cell_subset_locations(:,9)) median(cell_subset_locations(:,9)) mean(cell_subset_locations(:,10)) std(cell_subset_locations(:,10)) median(cell_subset_locations(:,10))];
        else
        header = {'Location Number' 'Z Slice' 'Local Pixel Row' 'Local Pixel Column' 'Global Parallel to Flow (um)' 'Global Perpindicular to Flow (um)' 'Diameter' [CD45marker ' value'] [CKmarker ' value'] [TUBmarker ' Value'] [ARmarker ' Value'] 'flag'};
        counts=[counts; cell_count_ofthistype mean(cell_subset_locations(:,7)) std(cell_subset_locations(:,7))  median(cell_subset_locations(:,7)) mean(cell_subset_locations(:,8)) std(cell_subset_locations(:,8)) median(cell_subset_locations(:,8)) mean(cell_subset_locations(:,9)) std(cell_subset_locations(:,9)) median(cell_subset_locations(:,9)) mean(cell_subset_locations(:,10)) std(cell_subset_locations(:,10)) median(cell_subset_locations(:,10)) mean(cell_subset_locations(:,11)) std(cell_subset_locations(:,11)) median(cell_subset_locations(:,11))];
        end;
     if(validation~=1)  
     xlswrite(strcat(dir,'CTC_Quant_locations_',day_for_filename,'_',time_for_filename,filename_corrected_tag), cell_count, namevec,'A1')
     else
     xlswrite(strcat(dir,'MANUAL_' ,xlsfilename), cell_count, namevec,'A1')
     end
    % stuff to use in the bar plot later
    list_of_namevecs(ii) = {namevec};
    list_of_cell_counts(ii) = cell_count_ofthistype;
  %cell array of 1 by 6
  if(validation~=1)  
    xlswrite(strcat(dir,'CTC_Quant_locations_',day_for_filename,'_',time_for_filename,filename_corrected_tag), header, namevec,'A2') % by defualt starts from A1
  else
     xlswrite(strcat(dir,'MANUAL_' , xlsfilename),  header, namevec,'A2')
  end
    if ~isempty(cell_subset_locations)
        if(validation~=1)
        xlswrite(strcat(dir,'CTC_Quant_locations_',day_for_filename,'_',time_for_filename,filename_corrected_tag), cell_subset_locations, namevec,'A3') % array under the header.
        else
     xlswrite(strcat(dir,'MANUAL_' , xlsfilename),  cell_subset_locations, namevec,'A3')
        end
   end
end
name_counts=[name_counts;list_of_namevecs];
