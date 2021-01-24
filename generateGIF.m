clear;
name='slope'; %gif�ļ�����
FPS=10; %ÿ��֡����
indexBegin=1; %��ʱ�ļ���ʼ��š�
indexEnd=400; %��ʱ�ļ�������š�
interval=10; %��ʱ�ļ������
frameNumber=numel(indexBegin:interval:indexEnd);
mov=moviein(frameNumber);
close(gcf);
tempFileName='temp/slope_temp'; %��ʱ�ļ����·������ֹ�����֮ǰ����
count=0;
for index=indexBegin:interval:indexEnd
    load([tempFileName,num2str(index),'.mat']); %����ʱ�ļ���
    view=View(model);
    view.plot('groupId'); %��ͼ����ʾ��Ԫ������š�
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
    disp([num2str(count),'/',num2str(frameNumber)]); %��ʾ���ȡ�
end