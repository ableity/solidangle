% writed by fanzhen meng
% 2019/10/09
% matlab R2018b
% The distance-dirven calculated by the SRM calculation for PET system
% IN 3D mode
% The LORs were coded in planar mode
% The plane YoZ is parallel to the horizontal plane
%                   D1
%                                 |Y
%                     --------- |----------->Z(
%                                 |
%                   D2
% The x axis is vertial to the horizontal plane
%
% voxel code : x ,y, z
% lor code: y z
% have finished and non-errors
clear all;
close all;
clc

%  *******************************************
% Detector parameter
CryNumY=77; CryNumZ=104;
%CryNumY=78; CryNumZ=104;
CrySize=[26 4 4]; %晶体尺寸？--------by李蕾
CryCoorY=-CryNumY*CrySize(2)/2+CrySize(2)/2:CrySize(2):CryNumY*CrySize(2)/2-CrySize(2)/2; % Y轴每个晶体的中心    --------  by李蕾
CryCoorZ=-CryNumZ*CrySize(2)/2+CrySize(2)/2:CrySize(2):CryNumZ*CrySize(2)/2-CrySize(2)/2; % Z轴每个晶体的中心    --------  by李蕾
CryNumPerHead=CryNumY*CryNumZ;%晶体总数？--------by李蕾


Dis=240;%mm the distance between the two detector heads
LORNum=CryNumY^2*CryNumZ^2;%LOR总数   --------  by李蕾
VoxNumPerCry=4;% the voxel number for one crystal
% Voxel Parameter
VoxNumX=240;
%VoxNumY=412; VoxNumZ=516;
VoxNumY=308; VoxNumZ=416;%四倍的晶体尺寸--------  by李蕾
VoxNumYZ=VoxNumY*VoxNumZ;%yz平面平行于平板的平面

VoxSize=1;%mm
VoxCoorX=-VoxNumX*VoxSize/2+VoxSize/2:VoxSize:VoxNumX*VoxSize/2-VoxSize/2;%计算体素中心   --------  by李蕾
VoxCoorY=-VoxNumY*VoxSize/2 +VoxSize/2:VoxSize:VoxNumY*VoxSize/2-VoxSize/2;
VoxCoorZ=-VoxNumZ*VoxSize/2+VoxSize/2:VoxSize:VoxNumZ*VoxSize/2-VoxSize/2;

gap=0.22;
VoxNum=VoxNumYZ*VoxNumX;

%G:/task8-DD/reconstruction/System matrix/data_80_four/SRM_DD
nWSaveName = sprintf('number_72_105crystal_atten_Delta3_Method4.raw');
fod_nW = fopen(nWSaveName,'wb');

voxelIndexSaveName = sprintf('indvoxel_72_105crystal_atten_Delta3_Method4.raw');
fod_voxelIndex = fopen(voxelIndexSaveName,'wb');

weightValueSaveName = sprintf('weight_72_105crystal_atten_Delta3_Method4.raw');
fod_weightValue = fopen(weightValueSaveName,'wb');

%% load MC simulation data 
%MC矩阵用来算doi   --------  by李蕾
load('./nonzero_ratio.mat');
nonzero_ratio=nonzero_ratio(13,:);
% Solid=zeros(1,VoxNum);
eps=10^-8;
theta=1;
%% anaylize parameters
% method I:  equal distance
DeltaCry=12;
% NumDiv=fix(CrySize(1)/DeltaCry);
% RatioPerDiv=fix(CrySize(1)/(2*NumDiv));
% DeltaWeight=zeros(1,NumDiv);
% for ind=1:NumDiv
%     DeltaWeight(ind)=sum(nonzero_ratio((ind-1)*RatioPerDiv+1:ind*RatioPerDiv));
% end
% DeltaWeight(end)=DeltaWeight(end)+nonzero_ratio(end);
% %DeepLen=2*[0 2 4 6 8 10];
% DeepLen=2*[0 3 6 9 ];
% % DeepLen=2*[0 6];
% method II: equal weight   
% DeltaWeight=[0.5 0.5];
% DeepLen=2*[0 5.33];
% DeltaWeight=[0.25 0.25 0.25 0.25];
% DeepLen=2*[0 2.7 5.33 8.57];

% DeltaWeight=[1/6 1/6 1/6 1/6 1/6 1/6];
% DeepLen=2*[0 1.9 3.6  5.44  7.57 10.17];

% 
% % method III:  unequal distance and unequal proportity 
% DeltaWeight=[0.28 0.72];
% DeepLen=2*[0  3]; 
% DeltaWeight=[0.18 0.82];
% DeepLen=2*[0  2]; 

DeltaWeight=[sum(nonzero_ratio(1))  sum(nonzero_ratio(2:3))  sum(nonzero_ratio(5:7))  sum(nonzero_ratio(8:13))];
DeepLen=2*[0 1 3 7];
% DeltaWeight=[0.18  0.20  0.34  0.32];
% DeepLen=2*[0 2 4 8];

% DeltaWeight=[0.08  0.1  0.2  0.18  0.22  0.22];
% DeepLen=2*[0 1 2 4 6 9];


OffAbandon=0;%abandon the slices in the end of the crystals
Start=0;

%%
norm=zeros(LORNum,1);
tic
u_LYSO=0.087;
coeff= 1;

for LORi=1%:CryNumZ
    for LORj=1%:CryNumY% the crystal in the up head
        
        for LORm=1:CryNumZ%blocki%((blocki-1)*(CryNumZ/blockNum)+1):blocki*CryNumZ/blockNum
            
            for LORn=1:CryNumY% the crystal in the down head
                
                [LORm LORn]
                P=zeros(1,VoxNum);
                tmp=zeros(1,VoxNum);
                Solid=zeros(1,VoxNum);
                LORIndex=LORn+(LORm-1)*CryNumY;% 下板LOR点的编号--------  by李蕾
                
                %% Depth in crystals
                for IndDoiUp=(Start+1):size(DeepLen,2)-OffAbandon  %1-4   --------  by李蕾
                        weightup=DeltaWeight(IndDoiUp);    %根据1-4的不同，找到MC中不同深度的DOI权重   --------  by李蕾
                        
                    for  IndDoiDown=(Start+1):size(DeepLen,2)-OffAbandon %双平板，up+down同上    --------  by李蕾
                         weightdown=DeltaWeight(IndDoiDown); %同上     --------  by李蕾
                  
                        OffsetUP=DeepLen(IndDoiUp);  %可能是深度值？    --------  by李蕾
                        OffsetDown=DeepLen(IndDoiDown);
                        
                        LORUp=[Dis/2+DeepLen(IndDoiUp) CryCoorY(LORj) CryCoorZ(LORi)]; %LOR的上坐标（晶体中心），因为由doi所以x轴也有值，上平板的晶体是固定（1，1）的  --------  by李蕾
                        LORDown=[-Dis/2-DeepLen(IndDoiDown) CryCoorY(LORn) CryCoorZ(LORm)];% center point LOR 的下左边   --------  by李蕾
                        
                        LORUpLD=[Dis/2+DeepLen(IndDoiUp) CryCoorY(LORj)-CrySize(2)/2+VoxSize/2 CryCoorZ(LORi)-CrySize(2)/2+VoxSize/2];  %上板第一个体素正上方的晶体位置    --------  by李蕾
                        LORDownLD=[-Dis/2-DeepLen(IndDoiDown) CryCoorY(LORn)-CrySize(2)/2+VoxSize/2  CryCoorZ(LORm)-CrySize(2)/2+VoxSize/2];% zuo xia jiao 
                        
                        kx=LORDown(1)-LORUp(1); ky=LORDown(2)-LORUp(2);kz=LORDown(3)-LORUp(3);
                        vecLOR=[kx ky kz]; % LOR向量    --------  by李蕾
                        vecY=[0 1 0];
                        vecZ=[0 0 1];
                        
                        lenLOR=sqrt(kx^2+ky^2+kz^2); %LOR长度    --------  by李蕾
                        
                        angleY=acos(abs(ky)/lenLOR);  %LOR与Y轴夹角
                        angleZ=acos(abs(kz)/lenLOR);  %LOR与Z轴夹角
                        
                        AttenLenUp=OffsetUP/(cos(atan(Dis/lenLOR)));%  这个在算什么？如果是晶体内部的衰减距离，那atan应该是acos  --------  by李蕾
                        AttenLenDown=OffsetDown/(cos(atan(Dis/lenLOR)));% 上面的问题无所谓，因为这个值没有调用
                        
                        AttenWeighUp=exp(-u_LYSO*AttenLenUp);% transform in the line
                        AttenWeighDown=exp(-u_LYSO*AttenLenDown);%transform in the line
                        
                        sliceEff=1;
                        if ky==0&&kz==0% 最角落的点  --------  by李蕾
                            
                            IndexY=ceil((LORUp(2)-(CryCoorY(1)-CrySize(2)/2))/CrySize(2));%the index of crystal  (-152-(-152-4/2))/4 = 1 ceil向上取整
                            IndexZ=ceil((LORUp(3)-(CryCoorZ(1)-CrySize(3)/2))/CrySize(3));%1
                            
                            IndexVoxY=(IndexY-1)*VoxNumPerCry+50+1;%51 为啥加50？可能是视野范围小一点   这个值应该是最下面一层的第一个体素 --------  by李蕾
                            IndexVoxZ=(IndexZ-1)*VoxNumPerCry+50+1;%51
                            
                            for tmpXi=1:VoxNumX
                                
                                IndexTemp= tmpXi+(IndexVoxY-1)*VoxNumX+ (IndexVoxZ-1)*VoxNumY*VoxNumX; %先x再y再z  算每个晶体对应的16个体素中的第一个体素的编号   --------  by李蕾
                                
                                for Voxi=1:VoxNumPerCry%Z
                                    for Voxj=1:VoxNumPerCry%Y
                                        Index=IndexTemp+(Voxj-1)*VoxNumX+(Voxi-1)*VoxNumY*VoxNumX;%P的序号， 根据indexTemp每个体素对应的最终编号
                                        
                                        point=[VoxCoorX(tmpXi)  VoxCoorY(IndexVoxY+Voxj-1)  VoxCoorZ(IndexVoxZ+Voxi-1)];%   当前遍历的体素坐标
                                        centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP);
                                        
                                        theta= SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP);
                                        
                                        P(Index)=VoxSize^2*theta*sliceEff;
                                        
                                    end
                                end
                            end
                            
                            
                        elseif ky~=0||kz~=0
                            
                            X=VoxCoorX;%体素中心   --------  by李蕾
                            
                            t=(X- LORUpLD(1))/kx;%应该是按照体素x轴位置对平板间距进行归一化后的距离
                            Y=ky*t+ LORUpLD(2);% ？？？    --------  by李蕾
                            Z=kz*t+ LORUpLD(3);
                            %
                            Y=Y(find(Y>VoxCoorY(1)-VoxSize/2&Y<VoxCoorY(end)+VoxSize/2)); %貌似在算所有边界内的体素    --------  by李蕾
                            Z=Z(find(Z>VoxCoorZ(1)-VoxSize/2&Z<VoxCoorZ(end)+VoxSize/2));
                            if ~isempty(Y)&&~isempty(Z)
                                YZIndex=intersect(find(Y>VoxCoorY(1)-VoxSize/2&Y<VoxCoorY(end)+VoxSize/2),find(Z>VoxCoorZ(1)-VoxSize/2&Z<VoxCoorZ(end)+VoxSize/2)); %为啥它又算一次?    --------  by李蕾
                                X=X(YZIndex);
                                Y=Y(YZIndex);
                                Z=Z(YZIndex);
                                
                                IndexInX=ceil((X-(VoxCoorX(1)-VoxSize/2))/VoxSize);
                                IndexInY=ceil((Y-(VoxCoorY(1)-VoxSize/2))/VoxSize);
                                IndexInZ=ceil((Z-(VoxCoorZ(1)-VoxSize/2))/VoxSize);
                                %
                                VarInY=VoxSize-mod(Y,VoxSize);
                                VarInZ=VoxSize-mod(Z,VoxSize);
                                
                                Index=IndexInX+(IndexInY-1)*VoxNumX+(IndexInZ-1)*VoxNumX*VoxNumY;
                                % Index=IndexInY+(IndexInZ-1)*VoxNumY+(IndexInX-1)*VoxNumYZ;
                                %*****************************************************************************
                                %                               CASE11
                                %*****************************************************************************
                                
                                for slicei=1:VoxNumX
                                    
                                    if VarInY(slicei)<1&&VarInZ(slicei)<1
                                        for tmpi=1:VoxNumPerCry+1%Z
                                            for tmpj=1:VoxNumPerCry+1% Y
                                                
                                                IndexTmp=Index(slicei)+(tmpj-1)*VoxNumX+(tmpi-1)*VoxNumX*VoxNumY;
                                                
                                                point=[VoxCoorX(slicei) VoxCoorY(IndexInY(slicei)+tmpi-1) VoxCoorZ(IndexInZ(slicei)+tmpj-1)];
                                                centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP);
                                                
                                                theta= SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP);
                                                
                                                if (tmpi~=1||tmpi~=VoxNumPerCry+1)&&(tmpj~=1||tmpj~=VoxNumPerCry+1)
                                                    
                                                    % P(IndexTmp)=VoxSize*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif  tmpi==1&&(tmpj~=1||tmpj~=VoxNumPerCry+1)
                                                    
                                                    % P(IndexTmp)=VarInZ(slicei)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif  tmpi==VoxNumPerCry+1&&(tmpj~=1||tmpj~=VoxNumPerCry+1)
                                                    
                                                    %  P(IndexTmp)=mod(Z(slicei),VoxSize)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif( tmpi~=1||tmpi~=VoxNumPerCry+1)&&tmpj==1
                                                    
                                                    % P(IndexTmp)=VarInY(slicei)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif (tmpi~=1||tmpi~=VoxNumPerCry+1)&&tmpj==VoxNumPerCry+1
                                                    
                                                    %P(IndexTmp)=mod(Y(slicei),VoxSize)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif tmpi==1&&tmpj==1
                                                    
                                                    %P(IndexTmp)=VarInZ(slicei)*VarInY(slicei)*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif tmpi==1&&tmpj==VoxNumPerCry+1
                                                    
                                                    % P(IndexTmp)=VarInZ(slicei)*mod(Y(slicei),VoxSize)*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif tmpi==5&&tmpj==1
                                                    
                                                    %  P(IndexTmp)=VarInY(slicei)*mod(Z(slicei),VoxSize)*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                elseif tmpi==5&&tmpj==VoxNumPerCry+1
                                                    
                                                    %P(IndexTmp)=mod(Y(slicei),VoxSize) *mod(Z(slicei),VoxSize)*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                    
                                                end
                                                
                                            end
                                        end
                                        
                                        %*****************************************************************************
                                        %                                                                                 CASE12
                                        %*****************************************************************************
                                        
                                    elseif   VarInY(slicei)<1&&VarInZ(slicei)==1
                                        for tmpi=1:VoxNumPerCry%Z
                                            for tmpj=1:VoxNumPerCry+1%Y
                                                
                                                IndexTmp=Index(slicei)+(tmpj-1)*VoxNumX+(tmpi-1)*VoxNumX*VoxNumY;
                                                point=[VoxCoorX(slicei) VoxCoorY(IndexInY(slicei)+tmpi-1) VoxCoorZ(IndexInZ(slicei)+tmpj-1)];
                                                centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP);
                                                theta= SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP);
                                                if tmpj~=1||tmpj~=VoxNumPerCry+1
                                                    
                                                    % P(IndexTmp)=VoxSize*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                elseif  tmpj==1
                                                    
                                                    %  P(IndexTmp)=VarInY(slicei)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                elseif tmpj==VoxNumPerCry+1
                                                    
                                                    % P(IndexTmp)=mod(Y(slicei),VoxSize)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                end
                                                
                                            end
                                        end
                                        
                                        %*****************************************************************************
                                        %                                                                                CASE13
                                        %*****************************************************************************
                                        
                                    elseif VarInY(slicei)==1&&VarInZ(slicei)<1
                                        for tmpi=1:VoxNumPerCry+1%Z
                                            for tmpj=1:VoxNumPerCry%Y
                                                
                                                IndexTmp=Index(slicei)+(tmpj-1)*VoxNumX+(tmpi-1)*VoxNumX*VoxNumY;
                                                point=[VoxCoorX(slicei) VoxCoorY(IndexInY(slicei)+tmpi-1) VoxCoorZ(IndexInZ(slicei)+tmpj-1)];
                                                centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP);
                                                theta= SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP);
                                                
                                                if tmpi~=1||tmpj~=VoxNumPerCry+1
                                                    
                                                    % P(IndexTmp)=VoxSize*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                elseif  tmpi==1
                                                    
                                                    % P(IndexTmp)=VarInZ(slicei)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                elseif tmpi==VoxNumPerCry+1
                                                    
                                                    %P(IndexTmp)=mod(Z(slicei),VoxSize)*VoxSize*theta*sliceEff;
                                                    P(IndexTmp)=theta;
                                                end
                                                
                                            end
                                        end
                                        %*****************************************************************************
                                        %                                                                                    CASE 14
                                        %*****************************************************************************
                                        
                                    elseif VarInY(slicei)==1&&VarInZ(slicei)==1
                                        
                                        for tmpi=1:VoxNumPerCry%Z
                                            for tmpj=1:VoxNumPerCry%Y
                                                
                                                IndexTmp=Index(slicei)+(tmpj-1)*VoxNumX+(tmpi-1)*VoxNumX*VoxNumY;
                                                point=[VoxCoorX(slicei) VoxCoorY(IndexInY(slicei)+tmpi-1) VoxCoorZ(IndexInZ(slicei)+tmpj-1)];
                                                centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP);
                                                theta=SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP);
                                                P(IndexTmp)=VoxSize^2*theta*sliceEff;
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                            end
                            
                        end
                       
                        
                       % coeff= coefficient(OffsetUP,OffsetDown,lenLOR,Dis,gap,CrySize,u_LYSO);
                        
                         tmp=P*coeff* weightup*weightdown+tmp;
%                         P=P*weightup*weightdown;
                        %tmp=P+tmp;
                        %tmp=P*WeightByEvents(IndDoiUp)*WeightByEvents(IndDoiDown)+tmp;% MC simulation
                        %tmp=P*(1-AttenWeighUp)*(1-AttenWeighDown); % supposed that the photon progagation is a straight line in the crystal t
                        %tmp=P*AttenWeighUp*AttenWeighDown; % supposed that the photon progagation is a straight line in the crystal t2
                    end
                end
                
                voxelIndex=find(tmp~=0);
                weightValue=tmp(voxelIndex);
                nW=size(voxelIndex,2);
                
                
                fwrite(fod_nW,nW, 'int32');
                fwrite(fod_voxelIndex, voxelIndex, 'int32');
                fwrite(fod_weightValue, weightValue, 'float');
                
                
            end
        end
    end
end

fclose(fod_nW);
fclose(fod_voxelIndex);
fclose(fod_weightValue);



toc








