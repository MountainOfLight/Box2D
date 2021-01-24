classdef Model<handle
    properties(Constant)
        DEFAULT_MATERIAL=[1200,30e6,15e6,0.01,0.4,1e4]; %默认材料，六个参数依次为：
        %密度、
        %法向刚度（单位值，且与单元半径成正比）、
        %切向刚度（单位值，且与单元半径成正比）、
        %断裂位移（单位值，且与单元半径的平方成正比）、
        %内摩擦系数、
        %粘聚力（单位值，且与单元半径的平方成正比），
        %单位均为国际标准单位。
        DEFAULT_BORDER_STIFFNESS_COEFFICIENT=30e6; %默认边界刚度系数。
        DEFAULT_DAMPING_COEFFICIENT=10000; %默认阻尼系数。
    end
    properties
        averageRadius; %平均半径。
        radiusRange; %半径变化范围，最小半径=平均半径-半径变化范围/2，最大半径=平均半径+半径变化范围/2。
        boxWidth; %箱子宽度。
        boxHeight; %箱子高度。
        leftBorder; %左边界位置。
        rightBorder; %右边界位置。
        bottomBorder; %下边界位置。
        topBorder; %上边界位置。
        material; %材料矩阵，每行为一种材料，第一行为默认材料，每行从左至右依次为：密度、法向刚度、切向刚度、断裂位移、内摩擦系数、粘聚力。
        groupId; %单元所属组号。
        materialId; %单元材料号。
        density; %单元密度。
        r; %单元半径。
        m; %单元质量。
        x; %单元X坐标。
        y; %单元Y坐标。
        isLockedX; %是否锁定单元X坐标（布尔值），默认为false。
        isLockedY; %是否锁定单元Y坐标（布尔值），默认为false。
        x0; %单元X坐标初始值。
        y0; %单元Y坐标初始值。
        vx; %单元X方向速度。
        vy; %单元Y方向速度。
        ax; %单元X方向加速度。
        ay; %单元Y方向加速度。
        kn; %单元法向刚度。
        ks; %单元切向刚度。
        xf; %单元断裂位移。
        miu; %单元内摩擦系数。
        c; %单元粘聚力。
        fnx; %单元受到的来自其他单元的X方向法向力合力。
        fny; %单元受到的来自其他单元的Y方向法向力合力。
        borderStiffnessCoefficient; %边界刚度系数，默认为DEFAULT_BORDER_STIFFNESS_COEFFICIENT。
        borderForceX; %单元受到的来自边界的X方向法向力。
        borderForceY; %单元受到的来自边界的Y方向法向力。
        isShearOn; %是否计算切向力（布尔值），默认为false。
        fsx; %单元受到的来自其他单元的X方向切向力合力。
        fsy; %单元受到的来自其他单元的Y方向切向力合力。
        isFrictionOn; %是否计算摩擦力（布尔值），默认为false。
        ffx; %单元受到的来自其他单元的X方向摩擦力合力。
        ffy; %单元受到的来自其他单元的Y方向摩擦力合力。
        borderFrictionCoefficient; %边界摩擦系数，默认为0。
        borderFrictionX; %单元受到的来自边界的X方向摩擦力。
        borderFrictionY; %单元受到的来自边界的Y方向摩擦力。
        g; %重力加速度。
        Gx; %单元受到的X方向体力。
        Gy; %单元受到的Y方向体力（即重力）。
        dampingCoefficient; %单元阻尼系数（单位值，且与单元半径成正比），默认为r*DEFAULT_DAMPING_COEFFICIENT。
        fdx; %单元受到的X方向阻尼力。
        fdy; %单元受到的Y方向阻尼力。
        maximumAccumulatedDisplacement; %最大累计位移，默认为0.1*averageRadius。
        accumulatedDisplacementX; %单元在X方向上的累计位移。
        accumulatedDisplacementY; %单元在Y方向上的累计位移。
        expansionRate; %检索邻居单元时的“膨胀率”，默认为1.1。
        neighbour; %邻居矩阵。
        connection; %连接矩阵。
        connectionStrengthCoefficient; %连接强度系数矩阵，和connection一一对应，默认均为1，大于1连接增强，小于1连接减弱。
        t; %时间。
        dt; %时间步，默认为1e-4s。
        customParameter; %自定义参数。
    end
    methods(Static)
        function [x,y]=transformCoordinate(x0,y0,theta)
            %坐标变换函数，在旋转单元、计算切向力和摩擦力时使用。
            x=x0.*cos(theta)+y0.*sin(theta);
            y=y0.*cos(theta)-x0.*sin(theta);
        end
    end
    methods
        function model=Model(averageRadius0,radiusRange0,boxWidth0,boxHeight0)
            %构造函数，averageRadius0为平均半径，radiusRange0为半径变化范围，boxWidth0为箱子宽度，boxHeight0为箱子高度。
            model.averageRadius=averageRadius0;
            model.radiusRange=radiusRange0;
            model.boxWidth=boxWidth0;
            model.boxHeight=boxHeight0;
            model.leftBorder=0;
            model.rightBorder=model.boxWidth;
            model.bottomBorder=0;
            model.topBorder=model.boxHeight;
            model.material=Model.DEFAULT_MATERIAL;
            gridSize=2*(model.averageRadius+model.radiusRange/2);
            gridSizeX=model.boxWidth/floor(model.boxWidth/gridSize);
            gridSizeY=model.boxHeight/floor(model.boxHeight/gridSize);
            [X,Y]=meshgrid(0+gridSizeX/2:gridSizeX:model.boxWidth-gridSizeX/2,0+gridSizeY/2:gridSizeY:model.boxHeight-gridSizeY/2);
            x=reshape(X,numel(X),1);
            y=reshape(Y,numel(Y),1);
            model.x=x;
            model.y=y;
            model.isLockedX=false(size(model.x));
            model.isLockedY=false(size(model.y));
            model.x0=model.x;
            model.y0=model.y;
            model.r=model.averageRadius+model.radiusRange*rand(mean([length(model.x),length(model.y)]),1)-model.radiusRange/2;
            model.groupId=ones(size(model.r));
            model.materialId=ones(size(model.r));
            model.density=model.material(1,1)*ones(size(model.r));
            model.m=(model.density).*(4/3*pi*(model.r).^3);
            model.vx=zeros(size(model.x));
            model.vy=zeros(size(model.y));
            model.ax=zeros(size(model.x));
            model.ay=zeros(size(model.y));
            model.kn=model.r*model.material(1,2);
            model.ks=model.r*model.material(1,3);
            model.xf=(model.r).^2*model.material(1,4);
            model.miu=model.material(1,5)*ones(size(model.r));
            model.c=(model.r).^2*model.material(1,6);
            model.fnx=zeros(size(model.x));
            model.fny=zeros(size(model.y));
            model.borderStiffnessCoefficient=Model.DEFAULT_BORDER_STIFFNESS_COEFFICIENT;
            model.borderForceX=zeros(size(model.x));
            model.borderForceY=zeros(size(model.y));
            model.isShearOn=false;
            model.fsx=zeros(size(model.x));
            model.fsy=zeros(size(model.y));
            model.isFrictionOn=false;
            model.ffx=zeros(size(model.x));
            model.ffy=zeros(size(model.y));
            model.borderFrictionCoefficient=0;
            model.borderFrictionX=zeros(size(model.x));
            model.borderFrictionY=zeros(size(model.y));
            model.g=-9.8;
            model.Gx=zeros(size(model.x));
            model.Gy=model.m*model.g;
            model.dampingCoefficient=model.r*Model.DEFAULT_DAMPING_COEFFICIENT;
            model.fdx=zeros(size(model.x));
            model.fdy=zeros(size(model.y));
            model.maximumAccumulatedDisplacement=0.1*model.averageRadius;
            model.accumulatedDisplacementX=zeros(size(model.x));
            model.accumulatedDisplacementY=zeros(size(model.y));
            model.expansionRate=1.1;
            model.searchNeighbour();
            model.t=0;
            model.dt=1e-4;
            model.customParameter=[];
            disp('已建立模型：');
            disp(['平均半径 ',num2str(model.averageRadius),' m，']);
            disp(['半径变化范围 ',num2str(model.averageRadius-model.radiusRange/2),' - ',num2str(model.averageRadius+model.radiusRange/2),' m，']);
            disp(['箱子宽度 ',num2str(model.boxWidth),' m，']);
            disp(['箱子高度 ',num2str(model.boxHeight),' m，']);
            disp(['单元总数 ',num2str(length(model.r)),' 个。']);
        end
        function resetModel(model)
            %重置模型。
            model.x0=model.x;
            model.y0=model.y;
            model.t=0;
            disp('已重置模型。');
        end
        function moveBorder(model,border,distance)
            %移动边界，border可取'left'、'right'、'bottom'或'top'，distance为移动距离。
            if distance==0
                return;
            end
            switch border
                case 'left'
                    model.leftBorder=model.leftBorder+distance;
                    if distance<0
                        disp(['已将左边界向左移动 ',num2str(-distance),' m。']);
                    else
                        disp(['已将左边界向右移动 ',num2str(distance),' m。']);
                    end
                case 'right'
                    model.rightBorder=model.rightBorder+distance;
                    if distance<0
                        disp(['已将右边界向左移动 ',num2str(-distance),' m。']);
                    else
                        disp(['已将右边界向右移动 ',num2str(distance),' m。']);
                    end
                case 'bottom'
                    model.bottomBorder=model.bottomBorder+distance;
                    if distance<0
                        disp(['已将下边界向下移动 ',num2str(-distance),' m。']);
                    else
                        disp(['已将下边界向上移动 ',num2str(distance),' m。']);
                    end
                case 'top'
                    model.topBorder=model.topBorder+distance;
                    if distance<0
                        disp(['已将上边界向下移动 ',num2str(-distance),' m。']);
                    else
                        disp(['已将上边界向上移动 ',num2str(distance),' m。']);
                    end
            end
            model.boxWidth=model.rightBorder-model.leftBorder;
            model.boxHeight=model.topBorder-model.bottomBorder;
        end
        function mID=addMaterial(model,material1)
            %添加材料，mID为新增材料号，material1为新增材料的行向量，其元素依次为：密度、法向刚度、切向刚度、断裂位移、内摩擦系数、粘聚力。
            model.material=[model.material;material1];
            mID=size(model.material,1);
            disp(['已添加材料，当前共 ',num2str(size(model.material,1)),' 种材料。']);
        end
        function value=getMaterial(model,mID,property)
            %获得材料，value为材料某种性质的值，mID为材料号，property可取'density'、'kn'、'ks'、'xf'、'miu'或'c'。
            switch property
                case 'density'
                    value=model.material(mID,1);
                    disp([num2str(mID),' 号材料的密度为 ',num2str(value),' kg/m^3。']);
                case 'kn'
                    value=model.material(mID,2);
                    disp([num2str(mID),' 号材料的法向刚度为 ',num2str(value),' N/m。']);
                case 'ks'
                    value=model.material(mID,3);
                    disp([num2str(mID),' 号材料的切向刚度为 ',num2str(value),' N/m。']);
                case 'xf'
                    value=model.material(mID,4);
                    disp([num2str(mID),' 号材料的断裂位移为 ',num2str(value),' m。']);
                case 'miu'
                    value=model.material(mID,5);
                    disp([num2str(mID),' 号材料的内摩擦系数为 ',num2str(value),'。']);
                case 'c'
                    value=model.material(mID,6);
                    disp([num2str(mID),' 号材料的粘聚力为 ',num2str(value),' N。']);
            end
        end
        function setMaterial(model,mID,property,value)
            %设置材料，mID为材料号（不得为1，即不得为默认材料），property可取'density'、'kn'、'ks'、'xf'、'miu'或'c'，value为材料某种性质的值。
            if mID==1
                disp('不得修改默认材料！');
                return;
            end
            switch property
                case 'density'
                    model.material(mID,1)=value;
                    disp(['已将 ',num2str(mID),' 号材料的密度修改为 ',num2str(value),' kg/m^3。']);
                case 'kn'
                    model.material(mID,2)=value;
                    disp(['已将 ',num2str(mID),' 号材料的法向刚度修改为 ',num2str(value),' N/m。']);
                case 'ks'
                    model.material(mID,3)=value;
                    disp(['已将 ',num2str(mID),' 号材料的切向刚度修改为 ',num2str(value),' N/m。']);
                case 'xf'
                    model.material(mID,4)=value;
                    disp(['已将 ',num2str(mID),' 号材料的断裂位移修改为 ',num2str(value),' m。']);
                case 'miu'
                    model.material(mID,5)=value;
                    disp(['已将 ',num2str(mID),' 号材料的内摩擦系数修改为 ',num2str(value),'。']);
                case 'c'
                    model.material(mID,6)=value;
                    disp(['已将 ',num2str(mID),' 号材料的粘聚力修改为 ',num2str(value),' N。']);
            end
            eID=find(model.materialId==mID);
            model.updateMicroParameters(eID);
        end
        function printAllMaterials(model)
            %打印所有材料。
            for i=1:size(model.material,1)
                disp(['[',num2str(i),']',...
                    '   density: ',num2str(model.material(i,1),'%.4e'),' kg/m^3',...
                    '   kn: ',num2str(model.material(i,2),'%.4e'),' N/m',...
                    '   ks: ',num2str(model.material(i,3),'%.4e'),' N/m',...
                    '   xf: ',num2str(model.material(i,4),'%.4e'),' m',...
                    '   miu: ',num2str(model.material(i,5)),...
                    '   c: ',num2str(model.material(i,6),'%.4e'),' N']);
            end
            disp('已打印所有材料。');
        end
        function removeMaterial(model,mID)
            %删除材料，mID为材料号（不得为1，即不得为默认材料）。
            if mID==1
                disp('不得删除默认材料！');
                return;
            end
            model.material(mID,:)=[];
            eID=find(model.materialId==mID);
            model.materialId(eID)=1;
            model.updateMicroParameters(eID);
            disp(['已删除材料，当前共 ',num2str(size(model.material,1)),' 种材料。']);
        end
        function number=getElementNumber(model)
            %获得单元总数，number为单元总数。
            number=length(model.r);
            disp(['单元总数 ',num2str(number),' 个。']);
        end
        function eID=addElement(model,x1,y1,r1,gID,mID)
            %添加单元，eID为新增单元编号的列向量，x1为新增单元X坐标的列向量，y1为新增单元Y坐标的列向量，r1为新增单元半径的列向量，gID为新增单元所属组号，mID为新增单元材料号。
            eIDbegin=length(model.r)+1;
            model.groupId=[model.groupId;gID*ones(size(r1))];
            model.materialId=[model.materialId;mID*ones(size(r1))];
            model.density=[model.density;model.material(mID,1)*ones(size(r1))];
            model.r=[model.r;r1];
            model.m=[model.m;(model.material(mID,1)*ones(size(r1))).*(4/3*pi*(r1).^3)];
            model.x=[model.x;x1];
            model.y=[model.y;y1];
            model.isLockedX=[model.isLockedX;false(size(x1))];
            model.isLockedY=[model.isLockedY;false(size(y1))];
            model.x0=model.x;
            model.y0=model.y;
            model.vx=[model.vx;zeros(size(x1))];
            model.vy=[model.vy;zeros(size(y1))];
            model.ax=[model.ax;zeros(size(x1))];
            model.ay=[model.ay;zeros(size(y1))];
            model.kn=[model.kn;model.material(mID,2)*ones(size(r1))];
            model.ks=[model.ks;model.material(mID,3)*ones(size(r1))];
            model.xf=[model.xf;model.material(mID,4)*ones(size(r1))];
            model.miu=[model.miu;model.material(mID,5)*ones(size(r1))];
            model.c=[model.c;model.material(mID,6)*ones(size(r1))];
            model.fnx=[model.fnx;zeros(size(x1))];
            model.fny=[model.fny;zeros(size(y1))];
            model.borderForceX=[model.borderForceX;zeros(size(x1))];
            model.borderForceY=[model.borderForceY;zeros(size(y1))];
            model.fsx=[model.fsx;zeros(size(x1))];
            model.fsy=[model.fsy;zeros(size(y1))];
            model.ffx=[model.ffx;zeros(size(x1))];
            model.ffy=[model.ffy;zeros(size(y1))];
            model.borderFrictionX=[model.borderFrictionX;zeros(size(x1))];
            model.borderFrictionY=[model.borderFrictionY;zeros(size(y1))];
            model.Gx=[model.Gx;zeros(size(x1))];
            model.Gy=[model.Gy;(model.material(mID,1)*ones(size(r1))).*(4/3*pi*(r1).^3)*model.g];
            model.dampingCoefficient=[model.dampingCoefficient;10000*ones(size(r1))];
            model.fdx=[model.fdx;zeros(size(x1))];
            model.fdy=[model.fdy;zeros(size(y1))];
            model.accumulatedDisplacementX=[model.accumulatedDisplacementX;zeros(size(x1))];
            model.accumulatedDisplacementY=[model.accumulatedDisplacementY;zeros(size(y1))];
            model.neighbour=[model.neighbour;zeros(size(r1,1),size(model.neighbour,2))];
            model.connection=[model.connection;zeros(size(r1,1),size(model.connection,2))];
            model.connectionStrengthCoefficient=[model.connectionStrengthCoefficient;zeros(size(r1,1),size(model.connection,2))];
            model.searchNeighbour();
            eIDend=length(model.r);
            eID=(eIDbegin:eIDend)';
            disp(['已添加单元，当前共 ',num2str(length(model.r)),' 个单元。']);
        end
        function moveElement(model,eID,dx,dy)
            %移动单元，eID为单元编号的列向量，dx为沿X轴正方向移动的距离，dy为沿Y轴正方向移动的距离。
            model.x(eID)=model.x(eID)+dx;
            model.y(eID)=model.y(eID)+dy;
            model.searchNeighbour();
            disp(['已将单元沿X轴正方向移动 ',num2str(dx),' m，沿Y轴正方向移动 ',num2str(dy),' m。']);
        end
        function rotateElement(model,eID,dtheta)
            %旋转单元，eID为单元编号的列向量，dtheta为沿逆时针方向旋转的角度。
            tempx=model.x(eID);
            tempy=model.y(eID);
            cx=mean(tempx);
            cy=mean(tempy);
            tempx=tempx-cx;
            tempy=tempy-cy;
            [tempx,tempy]=Model.transformCoordinate(tempx,tempy,dtheta/180*pi);
            tempx=tempx+cx;
            tempy=tempy+cy;
            model.x(eID)=tempx;
            model.y(eID)=tempy;
            model.searchNeighbour();
            disp(['已将单元沿逆时针方向旋转 ',num2str(dtheta),' °。']);
        end
        function setWallElement(model,eID)
            %设置墙单元，eID为单元编号的列向量。
            model.lockElement(eID,'X');
            model.lockElement(eID,'Y');
            disp('已设置墙单元。');
        end
        function lockElement(model,eID,direction)
            %锁定单元，eID为单元编号的列向量，direction可取'X'或'Y'。
            switch direction
                case 'X'
                    model.isLockedX(eID)=true;
                case 'Y'
                    model.isLockedY(eID)=true;
            end
            disp('已锁定单元。');
        end
        function unlockElement(model,eID,direction)
            %解锁单元，eID为单元编号的列向量，direction可取'X'或'Y'。
            switch direction
                case 'X'
                    model.isLockedX(eID)=false;
                case 'Y'
                    model.isLockedY(eID)=false;
            end
            disp('已解锁单元。');
        end
        function removeElement(model,eID)
            %删除单元，eID为单元编号的列向量。
            model.groupId(eID)=[];
            model.materialId(eID)=[];
            model.density(eID)=[];
            model.r(eID)=[];
            model.m(eID)=[];
            model.x(eID)=[];
            model.y(eID)=[];
            model.isLockedX(eID)=[];
            model.isLockedY(eID)=[];
            model.x0(eID)=[];
            model.y0(eID)=[];
            model.vx(eID)=[];
            model.vy(eID)=[];
            model.ax(eID)=[];
            model.ay(eID)=[];
            model.kn(eID)=[];
            model.ks(eID)=[];
            model.xf(eID)=[];
            model.miu(eID)=[];
            model.c(eID)=[];
            model.fnx(eID)=[];
            model.fny(eID)=[];
            model.borderForceX(eID)=[];
            model.borderForceY(eID)=[];
            model.fsx(eID)=[];
            model.fsy(eID)=[];
            model.ffx(eID)=[];
            model.ffy(eID)=[];
            model.borderFrictionX(eID)=[];
            model.borderFrictionY(eID)=[];
            model.Gx(eID)=[];
            model.Gy(eID)=[];
            model.dampingCoefficient(eID)=[];
            model.fdx(eID)=[];
            model.fdy(eID)=[];
            model.accumulatedDisplacementX(eID)=[];
            model.accumulatedDisplacementY(eID)=[];
            model.neighbour(eID,:)=[];
            model.connection(eID,:)=[];
            model.connectionStrengthCoefficient(eID,:)=[];
            model.searchNeighbour();
            disp(['已删除单元，当前共 ',num2str(length(model.r)),' 个单元。']);
        end
        function eID=getElementId(model,gID)
            %获得单元编号，eID为单元编号的列向量，gID为组号。
            eID=find(model.groupId==gID);
            disp('已获得单元编号。');
        end
        function eID=getElementIdInPolygon(model,borderx,bordery)
            %获得多边形内部的单元编号，eID为单元编号的列向量，borderx为多边形顶点的X坐标，bordery为多边形顶点的Y坐标。
            [in,on]=inpolygon(model.x,model.y,borderx,bordery);
            eID=find(in|on);
            disp('已获得多边形内部的单元编号。');
        end
        function layeringModel(model,layer)
            %将模型由下至上分层，layer为层面坐标的矩阵，其中：
            %第一行为第一层（最下层）底面的X坐标，
            %第二行为第一层（最下层）底面的Y坐标，
            %第三行为第一层（最下层）顶面和第二层（次下层）底面的X坐标，
            %第四行为第一层（最下层）顶面和第二层（次下层）底面的Y坐标，
            %第五行为第二层（次下层）顶面和第三层底面的X坐标，
            %第六行为第二层（次下层）顶面和第三层底面的Y坐标，
            %…………
            %以此类推。
            model.groupId(:)=0;
            for i=1:size(layer,1)/2-1
                [in,on]=inpolygon(model.x,model.y,[layer(2*i-1,:),fliplr(layer(2*i+1,:)),layer(2*i-1,1)],[layer(2*i,:),fliplr(layer(2*i+2,:)),layer(2*i,1)]);
                eID=find(in|on);
                model.setGroupId(eID,i);
            end
            model.removeElement(find(model.groupId==0));
            disp('已将模型由下至上分层。');
        end
        function setGroupId(model,eID,gID)
            %设置单元所属组号，eID为单元编号的列向量，gID为组号。
            model.groupId(eID)=gID;
            disp('已设置单元所属组号。');
        end
        function setMaterialId(model,eID,mID)
            %设置单元材料号，eID为单元编号的列向量，mID为材料号。
            model.materialId(eID)=mID;
            model.updateMicroParameters(eID);
            disp('已设置单元材料号。');
        end
        function updateMicroParameters(model,eID)
            %更新微观参数，eID为单元编号的列向量。
            model.density(eID)=model.material(model.materialId(eID),1);
            model.m(eID)=(model.density(eID)).*(4/3*pi*(model.r(eID)).^3);
            model.Gy(eID)=model.m(eID)*model.g;
            model.kn(eID)=model.material(model.materialId(eID),2);
            model.ks(eID)=model.material(model.materialId(eID),3);
            model.xf(eID)=model.material(model.materialId(eID),4);
            model.miu(eID)=model.material(model.materialId(eID),5);
            model.c(eID)=model.material(model.materialId(eID),6);
            disp('已更新微观参数。');
        end
        function setBorderStiffnessCoefficient(model,borderStiffnessCoefficient0)
            %设置边界刚度系数，borderStiffnessCoefficient0为边界刚度系数。
            model.borderStiffnessCoefficient=borderStiffnessCoefficient0;
            disp('已设置边界刚度系数。');
        end
        function setDampingCoefficient(model,eID,dampingCoefficient0)
            %设置阻尼系数，eID为单元编号的列向量，dampingCoefficient0为基准阻尼系数，实际阻尼系数与单元半径r成正比。
            model.dampingCoefficient(eID)=model.r(eID)*dampingCoefficient0;
            disp('已设置阻尼系数。');
        end
        function searchNeighbour(model)
            %邻居检索函数。
            disp('正在检索邻居…………');
            tempr=model.r*model.expansionRate;
            tempx=model.x;
            tempy=model.y;
            tempneighbour=model.neighbour;
            tempconnection=model.connection;
            tempconnectionStrengthCoefficient=model.connectionStrengthCoefficient;
            if ~isempty(tempneighbour)&&~isempty(tempconnection)&&~isempty(tempconnectionStrengthCoefficient)
                neighbour0=tempneighbour;
                tempneighbour=zeros(length(tempr),size(tempneighbour,2));
                connection0=tempconnection;
                connectionStrengthCoefficient0=tempconnectionStrengthCoefficient;
            else
                tempneighbour=zeros(length(tempr),8);
            end
            M=size(tempneighbour,1);
            for i=1:M
                neighbourFilterX=(tempx-tempr>=tempx(i)-tempr(i)&tempx-tempr<=tempx(i)+tempr(i))|(tempx+tempr>=tempx(i)-tempr(i)&tempx+tempr<=tempx(i)+tempr(i))|(tempx-tempr<tempx(i)-tempr(i)&tempx+tempr>tempx(i)+tempr(i));
                neighbourFilterY=(tempy-tempr>=tempy(i)-tempr(i)&tempy-tempr<=tempy(i)+tempr(i))|(tempy+tempr>=tempy(i)-tempr(i)&tempy+tempr<=tempy(i)+tempr(i))|(tempy-tempr<tempy(i)-tempr(i)&tempy+tempr>tempy(i)+tempr(i));
                neighbourId=(find(neighbourFilterX&neighbourFilterY))';
                neighbourId=neighbourId(neighbourId~=i);
                Ni=length(neighbourId);
                if Ni>size(tempneighbour,2)
                    tempneighbour=[tempneighbour,zeros(M,Ni-size(tempneighbour,2))];
                end
                tempneighbour(i,1:Ni)=neighbourId;
            end
            tempneighbour=tempneighbour(:,1:ceil(find(tempneighbour~=0,1,'last')/M));
            N=size(tempneighbour,2);
            tempconnection=zeros(M,N);
            tempconnectionStrengthCoefficient=zeros(M,N);
            if exist('neighbour0','var')&&exist('connection0','var')
                for i=1:M
                    for j=1:N
                        if tempneighbour(i,j)==0
                            break;
                        else
                            j0=find(neighbour0(i,:)==tempneighbour(i,j));
                            if ~isempty(j0)
                                tempconnection(i,j)=connection0(i,j0);
                                tempconnectionStrengthCoefficient(i,j)=connectionStrengthCoefficient0(i,j0);
                            end
                        end
                    end
                end
            end
            model.neighbour=tempneighbour;
            model.connection=tempconnection;
            model.connectionStrengthCoefficient=tempconnectionStrengthCoefficient;
        end
        function connectInternalConnection(model,gID,varargin)
            %连接组内连接，gID为组号，varargin为连接强度系数（可以省略，省略时默认为1）。
            eID=find(model.groupId==gID);
            for i=1:length(eID)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID(i),j)==eID)==1
                            model.connection(eID(i),j)=1;
                            if nargin==3
                                model.connectionStrengthCoefficient(eID(i),j)=varargin{1};
                            else
                                model.connectionStrengthCoefficient(eID(i),j)=1;
                            end
                        end
                    end
                end
            end
            disp('已连接组内连接。');
        end
        function connectExternalConnection(model,gID1,gID2,varargin)
            %连接组间连接，gID1、gID2为组号，varargin为连接强度系数（可以省略，省略时默认为1）。
            eID1=find(model.groupId==gID1);
            eID2=find(model.groupId==gID2);
            for i=1:length(eID1)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID1(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID1(i),j)==eID2)==1
                            model.connection(eID1(i),j)=1;
                            if nargin==4
                                model.connectionStrengthCoefficient(eID1(i),j)=varargin{1};
                            else
                                model.connectionStrengthCoefficient(eID1(i),j)=1;
                            end
                        end
                    end
                end
            end
            disp('已连接组间连接。');
        end
        function connectAllExternalConnection(model,gID,varargin)
            %连接所有组间连接，gID为组号，varargin为连接强度系数（可以省略，省略时默认为1）。
            eID=find(model.groupId==gID);
            for i=1:length(eID)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID(i),j)==eID)==0
                            model.connection(eID(i),j)=1;
                            if nargin==3
                                model.connectionStrengthCoefficient(eID(i),j)=varargin{1};
                            else
                                model.connectionStrengthCoefficient(eID(i),j)=1;
                            end
                        end
                    end
                end
            end
            disp('已连接所有组间连接。');
        end
        function connectAllConnection(model,varargin)
            %连接所有连接，varargin为连接强度系数（可以省略，省略时默认为1）。
            for i=1:size(model.neighbour,1)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(i,j)==0
                        break;
                    else
                        model.connection(i,j)=1;
                        if nargin==2
                            model.connectionStrengthCoefficient(i,j)=varargin{1};
                        else
                            model.connectionStrengthCoefficient(i,j)=1;
                        end
                    end
                end
            end
            disp('已连接所有连接。');
        end
        function breakInternalConnection(model,gID)
            %断开组内连接，gID为组号。
            eID=find(model.groupId==gID);
            for i=1:length(eID)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID(i),j)==eID)==1
                            model.connection(eID(i),j)=0;
                            model.connectionStrengthCoefficient(eID(i),j)=0;
                        end
                    end
                end
            end
            disp('已断开组内连接。');
        end
        function breakExternalConnection(model,gID1,gID2)
            %断开组间连接，gID1、gID2为组号。
            eID1=find(model.groupId==gID1);
            eID2=find(model.groupId==gID2);
            for i=1:length(eID1)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID1(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID1(i),j)==eID2)==1
                            model.connection(eID1(i),j)=0;
                            model.connectionStrengthCoefficient(eID1(i),j)=0;
                        end
                    end
                end
            end
            disp('已断开组间连接。');
        end
        function breakAllExternalConnection(model,gID)
            %断开所有组间连接，gID为组号。
            eID=find(model.groupId==gID);
            for i=1:length(eID)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(eID(i),j)==0
                        break;
                    else
                        if max(model.neighbour(eID(i),j)==eID)==0
                            model.connection(eID(i),j)=0;
                            model.connectionStrengthCoefficient(eID(i),j)=0;
                        end
                    end
                end
            end
            disp('已断开所有组间连接。');
        end
        function breakAllConnection(model)
            %断开所有连接。
            for i=1:size(model.neighbour,1)
                for j=1:size(model.neighbour,2)
                    if model.neighbour(i,j)==0
                        break;
                    else
                        model.connection(i,j)=0;
                        model.connectionStrengthCoefficient(i,j)=0;
                    end
                end
            end
            disp('已断开所有连接。');
        end
        function balanceModel(model,varargin)
            %平衡模型，varargin为可选输入参数：
            %当省略，即为balanceModel()时，平衡一次，
            %当有一个，即为balanceModel(n)时，平衡n次，
            %当有两个，即为balanceModel(m,n)时，平衡m*n次，且每平衡n次显示一次进度。
            if nargin==1
                tempr=model.r;
                tempx=model.x;
                tempy=model.y;
                tempisLockedX=model.isLockedX;
                tempisLockedY=model.isLockedY;
                tempx0=model.x0;
                tempy0=model.y0;
                tempvx=model.vx;
                tempvy=model.vy;
                tempkn=model.kn;
                tempks=model.ks;
                tempxf=model.xf;
                tempmiu=model.miu;
                tempc=model.c;
                tempisShearOn=model.isShearOn;
                tempisFrictionOn=model.isFrictionOn;
                tempneighbour=model.neighbour;
                tempconnection=model.connection;
                tempconnectionStrengthCoefficient=model.connectionStrengthCoefficient;
                [M,N]=size(tempneighbour);
                neighbourr=zeros(M,N);
                neighbourdx=zeros(M,N);
                neighbourdy=zeros(M,N);
                neighbourdx0=zeros(M,N);
                neighbourdy0=zeros(M,N);
                neighbourkn=zeros(M,N);
                neighbourks=zeros(M,N);
                neighbourxf=zeros(M,N);
                neighbourmiu=zeros(M,N);
                neighbourc=zeros(M,N);
                neighbourfnx=zeros(M,N);
                neighbourfny=zeros(M,N);
                neighbourfsx=zeros(M,N);
                neighbourfsy=zeros(M,N);
                neighbourffx=zeros(M,N);
                neighbourffy=zeros(M,N);
                for i=1:M
                    if tempisLockedX(i)&&tempisLockedY(i)
                        continue;
                    else
                        for j=1:N
                            if tempneighbour(i,j)==0
                                break;
                            else
                                neighbourr(i,j)=tempr(tempneighbour(i,j))+tempr(i);
                                neighbourdx(i,j)=tempx(tempneighbour(i,j))-tempx(i);
                                neighbourdy(i,j)=tempy(tempneighbour(i,j))-tempy(i);
                                neighbourdx0(i,j)=tempx0(tempneighbour(i,j))-tempx0(i);
                                neighbourdy0(i,j)=tempy0(tempneighbour(i,j))-tempy0(i);
                                neighbourkn(i,j)=tempkn(i)*tempkn(tempneighbour(i,j))/(tempkn(i)+tempkn(tempneighbour(i,j)));
                                neighbourks(i,j)=tempks(i)*tempks(tempneighbour(i,j))/(tempks(i)+tempks(tempneighbour(i,j)));
                                if tempconnection(i,j)~=0
                                    neighbourxf(i,j)=min([tempxf(i),tempxf(tempneighbour(i,j))])*tempconnectionStrengthCoefficient(i,j);
                                end
                                neighbourmiu(i,j)=max([tempmiu(i),tempmiu(tempneighbour(i,j))]);
                                if tempconnection(i,j)~=0
                                    neighbourc(i,j)=min([tempc(i),tempc(tempneighbour(i,j))])*tempconnectionStrengthCoefficient(i,j);
                                end
                                xn=sqrt(neighbourdx(i,j)^2+neighbourdy(i,j)^2)-neighbourr(i,j);
                                if xn<0||tempconnection(i,j)~=0
                                    neighbourfnx(i,j)=neighbourkn(i,j)*xn*neighbourdx(i,j)/sqrt(neighbourdx(i,j)^2+neighbourdy(i,j)^2);
                                    neighbourfny(i,j)=neighbourkn(i,j)*xn*neighbourdy(i,j)/sqrt(neighbourdx(i,j)^2+neighbourdy(i,j)^2);
                                end
                                if tempisShearOn&&tempconnection(i,j)~=0
                                    alpha=atan(neighbourdy0(i,j)/neighbourdx0(i,j));
                                    [centerx,centery]=Model.transformCoordinate(tempx(i),tempy(i),alpha);
                                    [neighbourx,neighboury]=Model.transformCoordinate(tempx(tempneighbour(i,j)),tempy(tempneighbour(i,j)),alpha);
                                    xs=neighboury-centery;
                                    [neighbourfsx(i,j),neighbourfsy(i,j)]=Model.transformCoordinate(0,neighbourks(i,j)*xs,-alpha);
                                end
                                if tempconnection(i,j)~=0&&(xn>neighbourxf(i,j)||sqrt(neighbourfsx(i,j)^2+neighbourfsy(i,j)^2)>neighbourmiu(i,j)*sqrt(neighbourfnx(i,j)^2+neighbourfny(i,j)^2)+neighbourc(i,j))
                                    tempconnection(i,j)=0;
                                end
                                if tempisFrictionOn&&xn<0
                                    beta=atan(neighbourdy(i,j)/neighbourdx(i,j));
                                    [centervx,centervy]=Model.transformCoordinate(tempvx(i),tempvy(i),beta);
                                    [neighbourvx,neighbourvy]=Model.transformCoordinate(tempvx(tempneighbour(i,j)),tempvy(tempneighbour(i,j)),beta);
                                    if centervy<neighbourvy
                                        [neighbourffx(i,j),neighbourffy(i,j)]=Model.transformCoordinate(0,neighbourmiu(i,j)*sqrt(neighbourfnx(i,j)^2+neighbourfny(i,j)^2),-beta);
                                    end
                                    if centervy>neighbourvy
                                        [neighbourffx(i,j),neighbourffy(i,j)]=Model.transformCoordinate(0,-neighbourmiu(i,j)*sqrt(neighbourfnx(i,j)^2+neighbourfny(i,j)^2),-beta);
                                    end
                                end
                            end
                        end
                    end
                end
                model.fnx=sum(neighbourfnx,2);
                model.fny=sum(neighbourfny,2);
                model.borderForceX=zeros(M,1);
                model.borderForceY=zeros(M,1);
                leftBorderFilter=tempx-tempr<model.leftBorder;
                rightBorderFilter=tempx+tempr>model.rightBorder;
                bottomBorderFilter=tempy-tempr<model.bottomBorder;
                topBorderFilter=tempy+tempr>model.topBorder;
                leftBorderForce=model.borderStiffnessCoefficient*abs(tempr-(tempx-model.leftBorder));
                rightBorderForce=-model.borderStiffnessCoefficient*abs(tempr-(model.rightBorder-tempx));
                bottomBorderForce=model.borderStiffnessCoefficient*abs(tempr-(tempy-model.bottomBorder));
                topBorderForce=-model.borderStiffnessCoefficient*abs(tempr-(model.topBorder-tempy));
                model.borderForceX(leftBorderFilter)=leftBorderForce(leftBorderFilter);
                model.borderForceX(rightBorderFilter)=rightBorderForce(rightBorderFilter);
                model.borderForceY(bottomBorderFilter)=bottomBorderForce(bottomBorderFilter);
                model.borderForceY(topBorderFilter)=topBorderForce(topBorderFilter);
                model.fsx=sum(neighbourfsx,2);
                model.fsy=sum(neighbourfsy,2);
                model.ffx=sum(neighbourffx,2);
                model.ffy=sum(neighbourffy,2);
                model.borderFrictionX=zeros(M,1);
                model.borderFrictionY=zeros(M,1);
                if model.borderFrictionCoefficient~=0
                    model.borderFrictionX=-sign(model.vx)*model.borderFrictionCoefficient.*abs(model.borderForceY);
                    model.borderFrictionY=-sign(model.vy)*model.borderFrictionCoefficient.*abs(model.borderForceX);
                end
                model.fdx=-model.dampingCoefficient.*model.vx;
                model.fdy=-model.dampingCoefficient.*model.vy;
                model.ax=(model.fnx+model.borderForceX+model.fsx+model.ffx+model.borderFrictionX+model.Gx+model.fdx)./model.m;
                model.ay=(model.fny+model.borderForceY+model.fsy+model.ffy+model.borderFrictionY+model.Gy+model.fdy)./model.m;
                model.ax(tempisLockedX)=0;
                model.ay(tempisLockedY)=0;
                dx=model.vx*model.dt+0.5*model.ax*model.dt^2;
                dy=model.vy*model.dt+0.5*model.ay*model.dt^2;
                dx(tempisLockedX)=0;
                dy(tempisLockedY)=0;
                model.x=tempx+dx;
                model.y=tempy+dy;
                model.vx=model.vx+model.ax*model.dt;
                model.vy=model.vy+model.ay*model.dt;
                model.vx(tempisLockedX)=0;
                model.vy(tempisLockedY)=0;
                model.accumulatedDisplacementX=model.accumulatedDisplacementX+dx;
                model.accumulatedDisplacementY=model.accumulatedDisplacementY+dy;
                if max(abs(model.accumulatedDisplacementX))>model.maximumAccumulatedDisplacement||max(abs(model.accumulatedDisplacementY))>model.maximumAccumulatedDisplacement
                    model.searchNeighbour();
                    model.accumulatedDisplacementX=zeros(M,1);
                    model.accumulatedDisplacementY=zeros(M,1);
                end
                model.t=model.t+model.dt;
            end
            if nargin==2
                for i=1:varargin{1}
                    model.balanceModel();
                end
            end
            if nargin==3
                for i=1:varargin{1}
                    model.balanceModel(varargin{2});
                    disp([num2str(i),'/',num2str(varargin{1})]);
                end
            end
        end
        function gravitySedimentation(model,majorCycle,minorCycle)
            %重力沉积，平衡majorCycle*minorCycle次，且每平衡minorCycle次显示一次进度。
            v0max=0.1*model.averageRadius;
            model.vx=2*v0max*rand(size(model.x))-v0max;
            model.vy=2*v0max*rand(size(model.y))-v0max;
            model.balanceModel(majorCycle,minorCycle);
            disp('已重力沉积。');
        end
    end
end