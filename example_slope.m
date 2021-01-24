clear;

model=Model(1,0.3,200,150); %建立平均半径1m，半径变化范围0.3m，箱子宽度200m，箱子高度150m的模型。
majorCycle=300;
minorCycle=1000;
model.gravitySedimentation(majorCycle,minorCycle); %重力沉积，平衡300*1000次，每平衡1000次显示一次进度。
view=View(model);
view.plot('r'); %绘图，显示单元半径。
save('model/slope_model1.mat','model'); %保存模型文件。
pause;

slopetxt=load('slope/slope.txt'); %打开层面文件slope.txt。
model.layeringModel(slopetxt); %将模型由下至上分层。
model.moveBorder('top',-50); %将上边界向下移动50m，使得模型高度变成100m。
soil1txt=load('material/soil1.txt'); %打开材料文件soil1.txt。
mID1=model.addMaterial(soil1txt); %将材料soil1添加到模型中。
model.setMaterialId(model.getElementId(1),mID1); %将组1单元的材料设置成soil1。
model.setMaterialId(model.getElementId(3),mID1); %将组3单元的材料设置成soil1。
soil2txt=load('material/soil2.txt'); %同上。
mID2=model.addMaterial(soil2txt);
model.setMaterialId(model.getElementId(2),mID2);
view.plot('materialId'); %绘图，显示单元材料号。
model.setWallElement(model.getElementId(1)); %将组1单元设置成墙单元。
save('model/slope_model2.mat','model'); %保存模型文件。
pause;

model.resetModel(); %重置模型。
majorCycle=400;
minorCycle=1000;
%平衡400*1000次，每平衡1000次显示一次进度。
for i=1:majorCycle
    model.balanceModel(minorCycle);
    save(['temp/slope_temp',num2str(i),'.mat'],'model'); %保存临时文件。
    disp([num2str(i),'/',num2str(majorCycle)]);
end
view.plot('displacementX'); %绘图，显示单元X方向位移。
save('model/slope_model3.mat','model'); %保存模型文件。