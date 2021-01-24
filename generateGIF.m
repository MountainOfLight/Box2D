clear;
name='slope'; %gif文件名。
FPS=10; %每秒帧数。
indexBegin=1; %临时文件起始序号。
indexEnd=400; %临时文件终了序号。
interval=10; %临时文件间隔。
frameNumber=numel(indexBegin:interval:indexEnd);
mov=moviein(frameNumber);
close(gcf);
tempFileName='temp/slope_temp'; %临时文件相对路径（截止到序号之前）。
count=0;
for index=indexBegin:interval:indexEnd
    load([tempFileName,num2str(index),'.mat']); %打开临时文件。
    view=View(model);
    view.plot('groupId'); %绘图，显示单元所属组号。
    count=count+1;
    mov(count)=getframe(view.f);
    im=frame2im(mov(count));
    [A,map]=rgb2ind(im,256);
    if index==indexBegin
        imwrite(A,map,['gif/',name,'.gif'],'gif', 'Loopcount',inf,'DelayTime',1/FPS);
    else
        imwrite(A,map,['gif/',name,'.gif'],'gif','WriteMode','append','DelayTime',1/FPS);
    end
    close(view.f);
    disp([num2str(count),'/',num2str(frameNumber)]); %显示进度。
end