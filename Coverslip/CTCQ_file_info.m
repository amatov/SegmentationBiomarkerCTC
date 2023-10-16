
for zstep=1:99,
if zstep<10,
        zindex(zstep,:)=strcat('0',num2str(zstep));
    else
        zindex(zstep,:)=num2str(zstep);
    end;
end;  

for position_step=1:999,% define text names of positions.  assume <1000 positions.
    if position_step<10,
        pindex(position_step,:)=strcat('00',num2str(position_step));
    elseif position_step<100,
        pindex(position_step,:)=strcat('0',num2str(position_step));
    else
        pindex(position_step,:)=num2str(position_step);
    end;
end;
found_file=0;
position_step=1000;




if location==1,%WCMC
while found_file==0, % look through directory to find out how many positions
    position_step=position_step-1;
    dapifile=strcat(dir,'\',dapistart,'p000',pindex(position_step,:),'t00000001z0',zindex(1,:),DAPIchannel,'.TIF');
   found_file=exist(dapifile, 'file');
end;
tot_positions=position_step;% total number of position files
found_file=0;
zstep=100;
while found_file==0, % look through directory to find out how many z stacks
    zstep=zstep-1;
    dapifile=strcat(dir,'\',dapistart,'p000',pindex(position_step,:),'t00000001z0',zindex(zstep,:),DAPIchannel,'.TIF');
   found_file=exist(dapifile, 'file');
end;
z_stack_depth=zstep; %set # of z slices based on existing files



elseif location==2,%CU-Ith
    
    
    while found_file==0, % look through directory to find out how many positions
       position_step=position_step-1;
       dapifile=strcat([dir  dapistart '_' DAPIchannelTEXT '_L' num2str(position_step) '_Sum.lsm']);
      found_file=exist(dapifile, 'file');
    end;
    tot_positions=position_step;% total number of position files

    % now get z stack depth
    temp=tiffread(dapifile);
    temp2=size(temp);
    z_stack_depth=temp2(2);
end;

