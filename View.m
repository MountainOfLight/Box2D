classdef View<handle
    properties
        model; %模型。
        boxVisible; %是否显示箱子（布尔值），默认为true。
        boxLineWidth; %箱子的线宽，默认为1。
        elementOutlineVisible; %是否显示单元轮廓（布尔值），默认为false。
        elementOutlineLineWidth; %单元轮廓的线宽，默认为1。
        plottingPrecision; %绘图精度，可取1-5之间的整数，取1时精度最低，取5时精度最高。
        connectionLineWidth; %连接线的线宽，默认为1。
        f; %当前窗口。
        ax; %当前坐标轴。
    end
    methods
        function view=View(model0)
            %构造函数，model0为模型。
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
            %绘图，property可取：
            %'groupId'，单元所属组号，
            %'materialId'，单元材料号，
            %'density'，单元密度，
            %'r'，单元半径，
            %'d'，单元直径，
            %'m'，单元质量，
            %'isLockedX'，%单元X坐标是否被锁定，
            %'isLockedY'，%单元Y坐标是否被锁定，
            %'displacementX'，单元X方向位移，
            %'displacementY'，单元Y方向位移，
            %'vx'，单元X方向速度，
            %'vy'，单元Y方向速度，
            %'ax'，单元X方向加速度，
            %'ay'，单元Y方向加速度，
            %'kn'，单元法向刚度，
            %'ks'，单元切向刚度，
            %'xf'，单元断裂位移，
            %'miu'，单元内摩擦系数，
            %'c'，单元粘聚力，
            %'sigmaX'，单元X方向法向力，
            %'sigmaY'，单元Y方向法向力，
            %'tauX'，单元X方向切向力，
            %'tauY'，单元Y方向切向力，
            %'frictionX'，单元X方向摩擦力，
            %'frictionY'，单元Y方向摩擦力，
            %'gravityX'，单元X方向体力，
            %'gravityY'，单元Y方向体力，
            %'dampingForceX'，单元X方向阻尼力，
            %'dampingForceY'，单元Y方向阻尼力，
            %'connection'，单元之间的连接，
            %自定义参数。
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