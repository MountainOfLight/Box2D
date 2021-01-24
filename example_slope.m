clear;

model=Model(1,0.3,200,150); %����ƽ���뾶1m���뾶�仯��Χ0.3m�����ӿ��200m�����Ӹ߶�150m��ģ�͡�
majorCycle=300;
minorCycle=1000;
model.gravitySedimentation(majorCycle,minorCycle); %����������ƽ��300*1000�Σ�ÿƽ��1000����ʾһ�ν��ȡ�
view=View(model);
view.plot('r'); %��ͼ����ʾ��Ԫ�뾶��
save('model/slope_model1.mat','model'); %����ģ���ļ���
pause;

slopetxt=load('slope/slope.txt'); %�򿪲����ļ�slope.txt��
model.layeringModel(slopetxt); %��ģ���������Ϸֲ㡣
model.moveBorder('top',-50); %���ϱ߽������ƶ�50m��ʹ��ģ�͸߶ȱ��100m��
soil1txt=load('material/soil1.txt'); %�򿪲����ļ�soil1.txt��
mID1=model.addMaterial(soil1txt); %������soil1��ӵ�ģ���С�
model.setMaterialId(model.getElementId(1),mID1); %����1��Ԫ�Ĳ������ó�soil1��
model.setMaterialId(model.getElementId(3),mID1); %����3��Ԫ�Ĳ������ó�soil1��
soil2txt=load('material/soil2.txt'); %ͬ�ϡ�
mID2=model.addMaterial(soil2txt);
model.setMaterialId(model.getElementId(2),mID2);
view.plot('materialId'); %��ͼ����ʾ��Ԫ���Ϻš�
model.setWallElement(model.getElementId(1)); %����1��Ԫ���ó�ǽ��Ԫ��
save('model/slope_model2.mat','model'); %����ģ���ļ���
pause;

model.resetModel(); %����ģ�͡�
majorCycle=400;
minorCycle=1000;
%ƽ��400*1000�Σ�ÿƽ��1000����ʾһ�ν��ȡ�
for i=1:majorCycle
    model.balanceModel(minorCycle);
    save(['temp/slope_temp',num2str(i),'.mat'],'model'); %������ʱ�ļ���
    disp([num2str(i),'/',num2str(majorCycle)]);
end
view.plot('displacementX'); %��ͼ����ʾ��ԪX����λ�ơ�
save('model/slope_model3.mat','model'); %����ģ���ļ���