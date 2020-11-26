%function theta= SolidAngle3D3(point,LORUp,kx,ky,kz,CrySize,lenLOR,angleY,angleZ,Dis)
function theta= SolidAngle3D5(centerPoint,LORUp,kx,ky,kz,CrySize,angleY, angleZ, Dis,lenLOR,OffsetUP)


tUp=(centerPoint(1)-Dis/2-OffsetUP)/kx;
YUp=centerPoint(2)-ky*tUp;
ZUp=centerPoint(3)-kz*tUp;


lenPtoUp=sqrt((centerPoint(1)-Dis/2-OffsetUP)^2+(centerPoint(2)-YUp)^2+(centerPoint(3)-ZUp)^2);


lenYSide=abs(LORUp(2)+CrySize(2)/2-YUp);%
lenZSide=abs(LORUp(3)+CrySize(3)/2-ZUp);%%%%%%%%%%%%%%%%%%%%%%

RY=lenYSide*sin(angleY);
RZ=lenZSide*sin(angleZ);

lenProY=CrySize(2)*sin(angleY);
lenProZ=CrySize(3)*sin(angleZ);


RY=min(RY,lenProY-RY);
RZ=min(RZ,lenProZ-RZ);


LY=lenPtoUp-(lenYSide)*cos(angleY);
LZ=lenPtoUp-(lenZSide)*cos(angleZ);
tmpY=(CrySize(2)-lenYSide)*cos(angleY);
tmpZ=(CrySize(3)-lenZSide)*cos(angleZ);
%CASE1
if LY>0&&LZ>0
    
    LY1=lenLOR-lenPtoUp-tmpY;
    LZ1=lenLOR-lenPtoUp-tmpZ;
    %case11  the voxel is  in the parallelogram region
    if  LY1>0&&LZ1>0
        thetaY=atan(RY/LY)+atan(RY/LY1);
        thetaZ=atan(RZ/LZ)+atan(RZ/LZ1);%the planar angle is calculated by one  side of crystal
        
        maxLY=max(LY,LY1);
        thetaY=min(atan(lenProY/maxLY),thetaY);%atan(lenProY/maxLY the planar angle is calculated by two  sides of crystal
        maxLZ=max(LZ,LZ1);
        thetaZ=min(atan(lenProZ/maxLZ),thetaZ);
     %case12   the voxel is  in the triangle  region
    elseif LY1<=0&&LZ1<=0
        thetaY=atan(lenProY/LY);
        thetaZ=atan(lenProZ/LZ);

    elseif  LY1>0&&LZ1<=0
        thetaY=atan(RY/LY)+atan(RY/LY1);
        thetaZ=atan(lenProZ/LZ);
        
        maxLY=max(LY,LY1);
        thetaY=min(atan(lenProY/maxLY),thetaY);
        
    elseif LY1<=0&&LZ1>0
        thetaY=atan(lenProY/LY);
        thetaZ=atan(RZ/LZ)+atan(RZ/LZ1);
        maxLZ=max(LZ,LZ1);
        thetaZ=min(atan(lenProZ/maxLZ),thetaZ);
        
    end
    
    %CASE 2
elseif LY<=0&&LZ<=0
    thetaY=atan(lenProY/(lenLOR-tmpY-lenPtoUp));
    thetaZ=atan(lenProZ/(lenLOR-tmpZ-lenPtoUp));
    %CASE 3
elseif LY>0&&LZ<=0
    LY1=lenLOR-lenPtoUp-tmpY;
    if  LY1>0
        thetaY=atan(RY/LY)+atan(RY/LY1);
        
        maxLY=max(LY,LY1);
        thetaY=min(atan(lenProY/maxLY),thetaY);
        
    else
        thetaY=atan(lenProY/LY);
    end
    thetaZ=atan(lenProZ/(lenLOR-tmpZ-lenPtoUp));
    %CASE 4
elseif LY<=0&&LZ>0
    LZ1=lenLOR-lenPtoUp-tmpZ;
    if  LZ1>0
        thetaZ=atan(RZ/LZ)+atan(RZ/LZ1);
         maxLZ=max(LZ,LZ1);
        thetaZ=min(atan(lenProZ/maxLZ),thetaZ);
    else
        thetaZ=atan(lenProZ/LZ);
    end
    thetaY=atan(lenProY/(lenLOR-tmpY-lenPtoUp));
end

% thetaY=atan(RY/lenPtoUp)+atan(RY/(lenLOR-lenPtoUp));
% thetaZ=atan(RZ/lenPtoUp)+atan(RZ/(lenLOR-lenPtoUp));




%  thetaY=RY/LY+RY/(lenLOR-LY);
%  thetaZ=RZ/LZ+RZ/(lenLOR-LZ);
%
%  thetay=CrySize(2)/LY;
%  thetaz=CrySize(3)/LZ;





theta1=thetaY*thetaZ/(2*pi);

sliceeff=1;%VoxSize/sin(angle1);

theta=theta1*sliceeff;

end

