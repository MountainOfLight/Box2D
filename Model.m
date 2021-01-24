classdef Model<handle
    properties(Constant)
        DEFAULT_MATERIAL=[1200,30e6,15e6,0.01,0.4,1e4]; %Ĭ�ϲ��ϣ�������������Ϊ��
        %�ܶȡ�
        %����նȣ���λֵ�����뵥Ԫ�뾶�����ȣ���
        %����նȣ���λֵ�����뵥Ԫ�뾶�����ȣ���
        %����λ�ƣ���λֵ�����뵥Ԫ�뾶��ƽ�������ȣ���
        %��Ħ��ϵ����
        %ճ��������λֵ�����뵥Ԫ�뾶��ƽ�������ȣ���
        %��λ��Ϊ���ʱ�׼��λ��
        DEFAULT_BORDER_STIFFNESS_COEFFICIENT=30e6; %Ĭ�ϱ߽�ն�ϵ����
        DEFAULT_DAMPING_COEFFICIENT=10000; %Ĭ������ϵ����
    end
    properties
        averageRadius; %ƽ���뾶��
        radiusRange; %�뾶�仯��Χ����С�뾶=ƽ���뾶-�뾶�仯��Χ/2�����뾶=ƽ���뾶+�뾶�仯��Χ/2��
        boxWidth; %���ӿ�ȡ�
        boxHeight; %���Ӹ߶ȡ�
        leftBorder; %��߽�λ�á�
        rightBorder; %�ұ߽�λ�á�
        bottomBorder; %�±߽�λ�á�
        topBorder; %�ϱ߽�λ�á�
        material; %���Ͼ���ÿ��Ϊһ�ֲ��ϣ���һ��ΪĬ�ϲ��ϣ�ÿ�д�����������Ϊ���ܶȡ�����նȡ�����նȡ�����λ�ơ���Ħ��ϵ����ճ������
        groupId; %��Ԫ������š�
        materialId; %��Ԫ���Ϻš�
        density; %��Ԫ�ܶȡ�
        r; %��Ԫ�뾶��
        m; %��Ԫ������
        x; %��ԪX���ꡣ
        y; %��ԪY���ꡣ
        isLockedX; %�Ƿ�������ԪX���꣨����ֵ����Ĭ��Ϊfalse��
        isLockedY; %�Ƿ�������ԪY���꣨����ֵ����Ĭ��Ϊfalse��
        x0; %��ԪX�����ʼֵ��
        y0; %��ԪY�����ʼֵ��
        vx; %��ԪX�����ٶȡ�
        vy; %��ԪY�����ٶȡ�
        ax; %��ԪX������ٶȡ�
        ay; %��ԪY������ٶȡ�
        kn; %��Ԫ����նȡ�
        ks; %��Ԫ����նȡ�
        xf; %��Ԫ����λ�ơ�
        miu; %��Ԫ��Ħ��ϵ����
        c; %��Ԫճ������
        fnx; %��Ԫ�ܵ�������������Ԫ��X��������������
        fny; %��Ԫ�ܵ�������������Ԫ��Y��������������
        borderStiffnessCoefficient; %�߽�ն�ϵ����Ĭ��ΪDEFAULT_BORDER_STIFFNESS_COEFFICIENT��
        borderForceX; %��Ԫ�ܵ������Ա߽��X����������
        borderForceY; %��Ԫ�ܵ������Ա߽��Y����������
        isShearOn; %�Ƿ����������������ֵ����Ĭ��Ϊfalse��
        fsx; %��Ԫ�ܵ�������������Ԫ��X����������������
        fsy; %��Ԫ�ܵ�������������Ԫ��Y����������������
        isFrictionOn; %�Ƿ����Ħ����������ֵ����Ĭ��Ϊfalse��
        ffx; %��Ԫ�ܵ�������������Ԫ��X����Ħ����������
        ffy; %��Ԫ�ܵ�������������Ԫ��Y����Ħ����������
        borderFrictionCoefficient; %�߽�Ħ��ϵ����Ĭ��Ϊ0��
        borderFrictionX; %��Ԫ�ܵ������Ա߽��X����Ħ������
        borderFrictionY; %��Ԫ�ܵ������Ա߽��Y����Ħ������
        g; %�������ٶȡ�
        Gx; %��Ԫ�ܵ���X����������
        Gy; %��Ԫ�ܵ���Y��������������������
        dampingCoefficient; %��Ԫ����ϵ������λֵ�����뵥Ԫ�뾶�����ȣ���Ĭ��Ϊr*DEFAULT_DAMPING_COEFFICIENT��
        fdx; %��Ԫ�ܵ���X������������
        fdy; %��Ԫ�ܵ���Y������������
        maximumAccumulatedDisplacement; %����ۼ�λ�ƣ�Ĭ��Ϊ0.1*averageRadius��
        accumulatedDisplacementX; %��Ԫ��X�����ϵ��ۼ�λ�ơ�
        accumulatedDisplacementY; %��Ԫ��Y�����ϵ��ۼ�λ�ơ�
        expansionRate; %�����ھӵ�Ԫʱ�ġ������ʡ���Ĭ��Ϊ1.1��
        neighbour; %�ھӾ���
        connection; %���Ӿ���
        connectionStrengthCoefficient; %����ǿ��ϵ�����󣬺�connectionһһ��Ӧ��Ĭ�Ͼ�Ϊ1������1������ǿ��С��1���Ӽ�����
        t; %ʱ�䡣
        dt; %ʱ�䲽��Ĭ��Ϊ1e-4s��
        customParameter; %�Զ��������
    end
    methods(Static)
        function [x,y]=transformCoordinate(x0,y0,theta)
            %����任����������ת��Ԫ��������������Ħ����ʱʹ�á�
            x=x0.*cos(theta)+y0.*sin(theta);
            y=y0.*cos(theta)-x0.*sin(theta);
        end
    end
    methods
        function model=Model(averageRadius0,radiusRange0,boxWidth0,boxHeight0)
            %���캯����averageRadius0Ϊƽ���뾶��radiusRange0Ϊ�뾶�仯��Χ��boxWidth0Ϊ���ӿ�ȣ�boxHeight0Ϊ���Ӹ߶ȡ�
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
            disp('�ѽ���ģ�ͣ�');
            disp(['ƽ���뾶 ',num2str(model.averageRadius),' m��']);
            disp(['�뾶�仯��Χ ',num2str(model.averageRadius-model.radiusRange/2),' - ',num2str(model.averageRadius+model.radiusRange/2),' m��']);
            disp(['���ӿ�� ',num2str(model.boxWidth),' m��']);
            disp(['���Ӹ߶� ',num2str(model.boxHeight),' m��']);
            disp(['��Ԫ���� ',num2str(length(model.r)),' ����']);
        end
        function resetModel(model)
            %����ģ�͡�
            model.x0=model.x;
            model.y0=model.y;
            model.t=0;
            disp('������ģ�͡�');
        end
        function moveBorder(model,border,distance)
            %�ƶ��߽磬border��ȡ'left'��'right'��'bottom'��'top'��distanceΪ�ƶ����롣
            if distance==0
                return;
            end
            switch border
                case 'left'
                    model.leftBorder=model.leftBorder+distance;
                    if distance<0
                        disp(['�ѽ���߽������ƶ� ',num2str(-distance),' m��']);
                    else
                        disp(['�ѽ���߽������ƶ� ',num2str(distance),' m��']);
                    end
                case 'right'
                    model.rightBorder=model.rightBorder+distance;
                    if distance<0
                        disp(['�ѽ��ұ߽������ƶ� ',num2str(-distance),' m��']);
                    else
                        disp(['�ѽ��ұ߽������ƶ� ',num2str(distance),' m��']);
                    end
                case 'bottom'
                    model.bottomBorder=model.bottomBorder+distance;
                    if distance<0
                        disp(['�ѽ��±߽������ƶ� ',num2str(-distance),' m��']);
                    else
                        disp(['�ѽ��±߽������ƶ� ',num2str(distance),' m��']);
                    end
                case 'top'
                    model.topBorder=model.topBorder+distance;
                    if distance<0
                        disp(['�ѽ��ϱ߽������ƶ� ',num2str(-distance),' m��']);
                    else
                        disp(['�ѽ��ϱ߽������ƶ� ',num2str(distance),' m��']);
                    end
            end
            model.boxWidth=model.rightBorder-model.leftBorder;
            model.boxHeight=model.topBorder-model.bottomBorder;
        end
        function mID=addMaterial(model,material1)
            %��Ӳ��ϣ�mIDΪ�������Ϻţ�material1Ϊ�������ϵ�����������Ԫ������Ϊ���ܶȡ�����նȡ�����նȡ�����λ�ơ���Ħ��ϵ����ճ������
            model.material=[model.material;material1];
            mID=size(model.material,1);
            disp(['����Ӳ��ϣ���ǰ�� ',num2str(size(model.material,1)),' �ֲ��ϡ�']);
        end
        function value=getMaterial(model,mID,property)
            %��ò��ϣ�valueΪ����ĳ�����ʵ�ֵ��mIDΪ���Ϻţ�property��ȡ'density'��'kn'��'ks'��'xf'��'miu'��'c'��
            switch property
                case 'density'
                    value=model.material(mID,1);
                    disp([num2str(mID),' �Ų��ϵ��ܶ�Ϊ ',num2str(value),' kg/m^3��']);
                case 'kn'
                    value=model.material(mID,2);
                    disp([num2str(mID),' �Ų��ϵķ���ն�Ϊ ',num2str(value),' N/m��']);
                case 'ks'
                    value=model.material(mID,3);
                    disp([num2str(mID),' �Ų��ϵ�����ն�Ϊ ',num2str(value),' N/m��']);
                case 'xf'
                    value=model.material(mID,4);
                    disp([num2str(mID),' �Ų��ϵĶ���λ��Ϊ ',num2str(value),' m��']);
                case 'miu'
                    value=model.material(mID,5);
                    disp([num2str(mID),' �Ų��ϵ���Ħ��ϵ��Ϊ ',num2str(value),'��']);
                case 'c'
                    value=model.material(mID,6);
                    disp([num2str(mID),' �Ų��ϵ�ճ����Ϊ ',num2str(value),' N��']);
            end
        end
        function setMaterial(model,mID,property,value)
            %���ò��ϣ�mIDΪ���Ϻţ�����Ϊ1��������ΪĬ�ϲ��ϣ���property��ȡ'density'��'kn'��'ks'��'xf'��'miu'��'c'��valueΪ����ĳ�����ʵ�ֵ��
            if mID==1
                disp('�����޸�Ĭ�ϲ��ϣ�');
                return;
            end
            switch property
                case 'density'
                    model.material(mID,1)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵ��ܶ��޸�Ϊ ',num2str(value),' kg/m^3��']);
                case 'kn'
                    model.material(mID,2)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵķ���ն��޸�Ϊ ',num2str(value),' N/m��']);
                case 'ks'
                    model.material(mID,3)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵ�����ն��޸�Ϊ ',num2str(value),' N/m��']);
                case 'xf'
                    model.material(mID,4)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵĶ���λ���޸�Ϊ ',num2str(value),' m��']);
                case 'miu'
                    model.material(mID,5)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵ���Ħ��ϵ���޸�Ϊ ',num2str(value),'��']);
                case 'c'
                    model.material(mID,6)=value;
                    disp(['�ѽ� ',num2str(mID),' �Ų��ϵ�ճ�����޸�Ϊ ',num2str(value),' N��']);
            end
            eID=find(model.materialId==mID);
            model.updateMicroParameters(eID);
        end
        function printAllMaterials(model)
            %��ӡ���в��ϡ�
            for i=1:size(model.material,1)
                disp(['[',num2str(i),']',...
                    '   density: ',num2str(model.material(i,1),'%.4e'),' kg/m^3',...
                    '   kn: ',num2str(model.material(i,2),'%.4e'),' N/m',...
                    '   ks: ',num2str(model.material(i,3),'%.4e'),' N/m',...
                    '   xf: ',num2str(model.material(i,4),'%.4e'),' m',...
                    '   miu: ',num2str(model.material(i,5)),...
                    '   c: ',num2str(model.material(i,6),'%.4e'),' N']);
            end
            disp('�Ѵ�ӡ���в��ϡ�');
        end
        function removeMaterial(model,mID)
            %ɾ�����ϣ�mIDΪ���Ϻţ�����Ϊ1��������ΪĬ�ϲ��ϣ���
            if mID==1
                disp('����ɾ��Ĭ�ϲ��ϣ�');
                return;
            end
            model.material(mID,:)=[];
            eID=find(model.materialId==mID);
            model.materialId(eID)=1;
            model.updateMicroParameters(eID);
            disp(['��ɾ�����ϣ���ǰ�� ',num2str(size(model.material,1)),' �ֲ��ϡ�']);
        end
        function number=getElementNumber(model)
            %��õ�Ԫ������numberΪ��Ԫ������
            number=length(model.r);
            disp(['��Ԫ���� ',num2str(number),' ����']);
        end
        function eID=addElement(model,x1,y1,r1,gID,mID)
            %��ӵ�Ԫ��eIDΪ������Ԫ��ŵ���������x1Ϊ������ԪX�������������y1Ϊ������ԪY�������������r1Ϊ������Ԫ�뾶����������gIDΪ������Ԫ������ţ�mIDΪ������Ԫ���Ϻš�
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
            disp(['����ӵ�Ԫ����ǰ�� ',num2str(length(model.r)),' ����Ԫ��']);
        end
        function moveElement(model,eID,dx,dy)
            %�ƶ���Ԫ��eIDΪ��Ԫ��ŵ���������dxΪ��X���������ƶ��ľ��룬dyΪ��Y���������ƶ��ľ��롣
            model.x(eID)=model.x(eID)+dx;
            model.y(eID)=model.y(eID)+dy;
            model.searchNeighbour();
            disp(['�ѽ���Ԫ��X���������ƶ� ',num2str(dx),' m����Y���������ƶ� ',num2str(dy),' m��']);
        end
        function rotateElement(model,eID,dtheta)
            %��ת��Ԫ��eIDΪ��Ԫ��ŵ���������dthetaΪ����ʱ�뷽����ת�ĽǶȡ�
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
            disp(['�ѽ���Ԫ����ʱ�뷽����ת ',num2str(dtheta),' �㡣']);
        end
        function setWallElement(model,eID)
            %����ǽ��Ԫ��eIDΪ��Ԫ��ŵ���������
            model.lockElement(eID,'X');
            model.lockElement(eID,'Y');
            disp('������ǽ��Ԫ��');
        end
        function lockElement(model,eID,direction)
            %������Ԫ��eIDΪ��Ԫ��ŵ���������direction��ȡ'X'��'Y'��
            switch direction
                case 'X'
                    model.isLockedX(eID)=true;
                case 'Y'
                    model.isLockedY(eID)=true;
            end
            disp('��������Ԫ��');
        end
        function unlockElement(model,eID,direction)
            %������Ԫ��eIDΪ��Ԫ��ŵ���������direction��ȡ'X'��'Y'��
            switch direction
                case 'X'
                    model.isLockedX(eID)=false;
                case 'Y'
                    model.isLockedY(eID)=false;
            end
            disp('�ѽ�����Ԫ��');
        end
        function removeElement(model,eID)
            %ɾ����Ԫ��eIDΪ��Ԫ��ŵ���������
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
            disp(['��ɾ����Ԫ����ǰ�� ',num2str(length(model.r)),' ����Ԫ��']);
        end
        function eID=getElementId(model,gID)
            %��õ�Ԫ��ţ�eIDΪ��Ԫ��ŵ���������gIDΪ��š�
            eID=find(model.groupId==gID);
            disp('�ѻ�õ�Ԫ��š�');
        end
        function eID=getElementIdInPolygon(model,borderx,bordery)
            %��ö�����ڲ��ĵ�Ԫ��ţ�eIDΪ��Ԫ��ŵ���������borderxΪ����ζ����X���꣬borderyΪ����ζ����Y���ꡣ
            [in,on]=inpolygon(model.x,model.y,borderx,bordery);
            eID=find(in|on);
            disp('�ѻ�ö�����ڲ��ĵ�Ԫ��š�');
        end
        function layeringModel(model,layer)
            %��ģ���������Ϸֲ㣬layerΪ��������ľ������У�
            %��һ��Ϊ��һ�㣨���²㣩�����X���꣬
            %�ڶ���Ϊ��һ�㣨���²㣩�����Y���꣬
            %������Ϊ��һ�㣨���²㣩����͵ڶ��㣨���²㣩�����X���꣬
            %������Ϊ��һ�㣨���²㣩����͵ڶ��㣨���²㣩�����Y���꣬
            %������Ϊ�ڶ��㣨���²㣩����͵���������X���꣬
            %������Ϊ�ڶ��㣨���²㣩����͵���������Y���꣬
            %��������
            %�Դ����ơ�
            model.groupId(:)=0;
            for i=1:size(layer,1)/2-1
                [in,on]=inpolygon(model.x,model.y,[layer(2*i-1,:),fliplr(layer(2*i+1,:)),layer(2*i-1,1)],[layer(2*i,:),fliplr(layer(2*i+2,:)),layer(2*i,1)]);
                eID=find(in|on);
                model.setGroupId(eID,i);
            end
            model.removeElement(find(model.groupId==0));
            disp('�ѽ�ģ���������Ϸֲ㡣');
        end
        function setGroupId(model,eID,gID)
            %���õ�Ԫ������ţ�eIDΪ��Ԫ��ŵ���������gIDΪ��š�
            model.groupId(eID)=gID;
            disp('�����õ�Ԫ������š�');
        end
        function setMaterialId(model,eID,mID)
            %���õ�Ԫ���Ϻţ�eIDΪ��Ԫ��ŵ���������mIDΪ���Ϻš�
            model.materialId(eID)=mID;
            model.updateMicroParameters(eID);
            disp('�����õ�Ԫ���Ϻš�');
        end
        function updateMicroParameters(model,eID)
            %����΢�۲�����eIDΪ��Ԫ��ŵ���������
            model.density(eID)=model.material(model.materialId(eID),1);
            model.m(eID)=(model.density(eID)).*(4/3*pi*(model.r(eID)).^3);
            model.Gy(eID)=model.m(eID)*model.g;
            model.kn(eID)=model.material(model.materialId(eID),2);
            model.ks(eID)=model.material(model.materialId(eID),3);
            model.xf(eID)=model.material(model.materialId(eID),4);
            model.miu(eID)=model.material(model.materialId(eID),5);
            model.c(eID)=model.material(model.materialId(eID),6);
            disp('�Ѹ���΢�۲�����');
        end
        function setBorderStiffnessCoefficient(model,borderStiffnessCoefficient0)
            %���ñ߽�ն�ϵ����borderStiffnessCoefficient0Ϊ�߽�ն�ϵ����
            model.borderStiffnessCoefficient=borderStiffnessCoefficient0;
            disp('�����ñ߽�ն�ϵ����');
        end
        function setDampingCoefficient(model,eID,dampingCoefficient0)
            %��������ϵ����eIDΪ��Ԫ��ŵ���������dampingCoefficient0Ϊ��׼����ϵ����ʵ������ϵ���뵥Ԫ�뾶r�����ȡ�
            model.dampingCoefficient(eID)=model.r(eID)*dampingCoefficient0;
            disp('����������ϵ����');
        end
        function searchNeighbour(model)
            %�ھӼ���������
            disp('���ڼ����ھӡ�������');
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
            %�����������ӣ�gIDΪ��ţ�vararginΪ����ǿ��ϵ��������ʡ�ԣ�ʡ��ʱĬ��Ϊ1����
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
            disp('�������������ӡ�');
        end
        function connectExternalConnection(model,gID1,gID2,varargin)
            %����������ӣ�gID1��gID2Ϊ��ţ�vararginΪ����ǿ��ϵ��������ʡ�ԣ�ʡ��ʱĬ��Ϊ1����
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
            disp('������������ӡ�');
        end
        function connectAllExternalConnection(model,gID,varargin)
            %��������������ӣ�gIDΪ��ţ�vararginΪ����ǿ��ϵ��������ʡ�ԣ�ʡ��ʱĬ��Ϊ1����
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
            disp('����������������ӡ�');
        end
        function connectAllConnection(model,varargin)
            %�����������ӣ�vararginΪ����ǿ��ϵ��������ʡ�ԣ�ʡ��ʱĬ��Ϊ1����
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
            disp('�������������ӡ�');
        end
        function breakInternalConnection(model,gID)
            %�Ͽ��������ӣ�gIDΪ��š�
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
            disp('�ѶϿ��������ӡ�');
        end
        function breakExternalConnection(model,gID1,gID2)
            %�Ͽ�������ӣ�gID1��gID2Ϊ��š�
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
            disp('�ѶϿ�������ӡ�');
        end
        function breakAllExternalConnection(model,gID)
            %�Ͽ�����������ӣ�gIDΪ��š�
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
            disp('�ѶϿ�����������ӡ�');
        end
        function breakAllConnection(model)
            %�Ͽ��������ӡ�
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
            disp('�ѶϿ��������ӡ�');
        end
        function balanceModel(model,varargin)
            %ƽ��ģ�ͣ�vararginΪ��ѡ���������
            %��ʡ�ԣ���ΪbalanceModel()ʱ��ƽ��һ�Σ�
            %����һ������ΪbalanceModel(n)ʱ��ƽ��n�Σ�
            %������������ΪbalanceModel(m,n)ʱ��ƽ��m*n�Σ���ÿƽ��n����ʾһ�ν��ȡ�
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
            %����������ƽ��majorCycle*minorCycle�Σ���ÿƽ��minorCycle����ʾһ�ν��ȡ�
            v0max=0.1*model.averageRadius;
            model.vx=2*v0max*rand(size(model.x))-v0max;
            model.vy=2*v0max*rand(size(model.y))-v0max;
            model.balanceModel(majorCycle,minorCycle);
            disp('������������');
        end
    end
end