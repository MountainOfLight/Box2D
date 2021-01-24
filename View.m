classdef View<handle
    properties
        model; %ģ�͡�
        boxVisible; %�Ƿ���ʾ���ӣ�����ֵ����Ĭ��Ϊtrue��
        boxLineWidth; %���ӵ��߿�Ĭ��Ϊ1��
        elementOutlineVisible; %�Ƿ���ʾ��Ԫ����������ֵ����Ĭ��Ϊfalse��
        elementOutlineLineWidth; %��Ԫ�������߿�Ĭ��Ϊ1��
        plottingPrecision; %��ͼ���ȣ���ȡ1-5֮���������ȡ1ʱ������ͣ�ȡ5ʱ������ߡ�
        connectionLineWidth; %�����ߵ��߿�Ĭ��Ϊ1��
        f; %��ǰ���ڡ�
        ax; %��ǰ�����ᡣ
    end
    methods
        function view=View(model0)
            %���캯����model0Ϊģ�͡�
            view.model=model0;
            view.boxVisible=true;
            view.boxLineWidth=1;
            view.elementOutlineVisible=false;
            view.elementOutlineLineWidth=1;
            view.plottingPrecision=3;
            view.connectionLineWidth=1;
            view.f=[];
            view.ax=[];
        end
        function plot(view,property)
            %��ͼ��property��ȡ��
            %'groupId'����Ԫ������ţ�
            %'materialId'����Ԫ���Ϻţ�
            %'density'����Ԫ�ܶȣ�
            %'r'����Ԫ�뾶��
            %'d'����Ԫֱ����
            %'m'����Ԫ������
            %'isLockedX'��%��ԪX�����Ƿ�������
            %'isLockedY'��%��ԪY�����Ƿ�������
            %'displacementX'����ԪX����λ�ƣ�
            %'displacementY'����ԪY����λ�ƣ�
            %'vx'����ԪX�����ٶȣ�
            %'vy'����ԪY�����ٶȣ�
            %'ax'����ԪX������ٶȣ�
            %'ay'����ԪY������ٶȣ�
            %'kn'����Ԫ����նȣ�
            %'ks'����Ԫ����նȣ�
            %'xf'����Ԫ����λ�ƣ�
            %'miu'����Ԫ��Ħ��ϵ����
            %'c'����Ԫճ������
            %'sigmaX'����ԪX����������
            %'sigmaY'����ԪY����������
            %'tauX'����ԪX������������
            %'tauY'����ԪY������������
            %'frictionX'����ԪX����Ħ������
            %'frictionY'����ԪY����Ħ������
            %'gravityX'����ԪX����������
            %'gravityY'����ԪY����������
            %'dampingForceX'����ԪX������������
            %'dampingForceY'����ԪY������������
            %'connection'����Ԫ֮������ӣ�
            %�Զ��������
            view.f=figure();
            set(view.f,'Color','w');
            view.ax=axes();
            axis(view.ax,'equal');
            hold(view.ax,'on');
            title(view.ax,property);
            xlabel(view.ax,'x(m)');
            ylabel(view.ax,'y(m)');
            colormap(view.ax,'jet');
            colorbar(view.ax,'eastoutside');
            theta=0:6*(6-view.plottingPrecision):360;
            x=view.model.x+view.model.r*cosd(theta);
            y=view.model.y+view.model.r*sind(theta);
            if view.boxVisible
                rectangle(view.ax,'Position',[view.model.leftBorder,view.model.bottomBorder,view.model.boxWidth,view.model.boxHeight],'LineStyle','-','LineWidth',view.boxLineWidth,'EdgeColor','k');
                axis(view.ax,[view.model.leftBorder,view.model.rightBorder,view.model.bottomBorder,view.model.topBorder]);
            else
                axis(view.ax,[min(x(:)),max(x(:)),min(y(:)),max(y(:))]);
            end
            switch property
                case 'groupId'
                    c=view.model.groupId;
                case 'materialId'
                    c=view.model.materialId;
                case 'density'
                    c=view.model.density;
                case 'r'
                    c=view.model.r;
                case 'd'
                    c=2*view.model.r;
                case 'm'
                    c=view.model.m;
                case 'isLockedX'
                    c=zeros(size(view.model.isLockedX));
                    c(view.model.isLockedX)=1;
                case 'isLockedY'
                    c=zeros(size(view.model.isLockedY));
                    c(view.model.isLockedY)=1;
                case 'displacementX'
                    c=view.model.x-view.model.x0;
                case 'displacementY'
                    c=view.model.y-view.model.y0;
                case 'vx'
                    c=view.model.vx;
                case 'vy'
                    c=view.model.vy;
                case 'ax'
                    c=view.model.ax;
                case 'ay'
                    c=view.model.ay;
                case 'kn'
                    c=view.model.kn;
                case 'ks'
                    c=view.model.ks;
                case 'xf'
                    c=view.model.xf;
                case 'miu'
                    c=view.model.miu;
                case 'c'
                    c=view.model.c;
                case 'sigmaX'
                    c=view.model.fnx+view.model.borderForceX;
                case 'sigmaY'
                    c=view.model.fny+view.model.borderForceY;
                case 'tauX'
                    c=view.model.fsx;
                case 'tauY'
                    c=view.model.fsy;
                case 'frictionX'
                    c=view.model.ffx+view.model.borderFrictionX;
                case 'frictionY'
                    c=view.model.ffy+view.model.borderFrictionY;
                case 'gravityX'
                    c=view.model.Gx;
                case 'gravityY'
                    c=view.model.Gy;
                case 'dampingForceX'
                    c=view.model.fdx;
                case 'dampingForceY'
                    c=view.model.fdy;
                otherwise
                    if isfield(view.model.customParameter,property)
                        c=getfield(view.model.customParameter,property);
                    end
            end
            if exist('c','var')
                for i=1:length(view.model.r)
                    if view.elementOutlineVisible
                        fill(view.ax,x(i,:),y(i,:),c(i),'LineStyle','-','LineWidth',view.elementOutlineLineWidth,'EdgeColor','k');
                    else
                        fill(view.ax,x(i,:),y(i,:),c(i),'LineStyle','none');
                    end
                end
            else
                if strcmp(property,'connection')
                    if view.elementOutlineVisible
                        for i=1:length(view.model.r)
                            plot(view.ax,x(i,:),y(i,:),'LineStyle','-','LineWidth',view.elementOutlineLineWidth,'Color','k');
                        end
                    end
                    isPlotted=true(size(view.model.connection));
                    isPlotted(logical(view.model.connection))=false;
                    for i=1:size(view.model.neighbour,1)
                        for j=1:size(view.model.neighbour,2)
                            if ~isPlotted(i,j)
                                if view.model.connectionStrengthCoefficient(i,j)<1
                                    connectionLineColor='r';
                                elseif view.model.connectionStrengthCoefficient(i,j)>1
                                    connectionLineColor='b';
                                else
                                    connectionLineColor='k';
                                end
                                plot(view.ax,[view.model.x(i),view.model.x(view.model.neighbour(i,j))],[view.model.y(i),view.model.y(view.model.neighbour(i,j))],'LineStyle','-','LineWidth',view.connectionLineWidth,'Color',connectionLineColor);
                                j0=find(view.model.neighbour(view.model.neighbour(i,j),:)==i);
                                isPlotted(view.model.neighbour(i,j),j0)=true;
                            end
                        end
                    end
                    colorbar(view.ax,'off');
                end
            end
        end
    end
end