function centerPoint= findCen(point,LORUp,kx,ky,kz,CrySize,VoxSize,Dis,OffsetUP)
centerPoint=point;

YSide=[LORUp(2)+CrySize(2)/2 LORUp(2)-CrySize(2)/2];
ZSide=[LORUp(3)+CrySize(3)/2 LORUp(3)-CrySize(3)/2];


tUp=(point(1)-Dis/2-OffsetUP)/kx;
YUp=point(2)-ky*tUp; 
ZUp=point(3)-kz*tUp;

if YUp>YSide(2)&&YUp<YSide(1)&&ZUp>ZSide(2)&&ZUp<ZSide(1)%
    centerPoint=point;
elseif YUp>YSide(2)&&YUp<YSide(1)&&ZUp>=ZSide(1)%
    Ztmp=kz*tUp+ZSide(1);
    centerPoint(3)=point(3)-0.5*VoxSize+0.5*(Ztmp-(point(3)-0.5*VoxSize));
elseif YUp>YSide(2)&&YUp<YSide(1)&&ZUp<=ZSide(2)%
    Ztmp=kz*tUp+ZSide(2);
    centerPoint(3)=point(3)+0.5*VoxSize-0.5*((point(3)+0.5*VoxSize)-Ztmp);
elseif ZUp>ZSide(2)&&ZUp<ZSide(1)&&YUp>=YSide(1)%
    Ytmp=ky*tUp+YSide(1);
    centerPoint(2)=point(2)-0.5*VoxSize+0.5*(Ytmp-(point(2)-0.5*VoxSize));
elseif ZUp>ZSide(2)&&ZUp<ZSide(1)&&YUp<=YSide(2)%
    Ytmp=ky*tUp+YSide(2);
    centerPoint(2)=point(2)+0.5*VoxSize-0.5*((point(2)+0.5*VoxSize)-Ytmp);
    
 elseif ZUp>=ZSide(1)&&YUp>=YSide(1)%
    Ytmp=ky*tUp+YSide(1);
    Ztmp=kz*tUp+ZSide(1);
    centerPoint(2)=point(2)-0.5*VoxSize+0.5*(Ytmp-(point(2)-0.5*VoxSize));
    centerPoint(3)=point(3)-0.5*VoxSize+0.5*(Ztmp-(point(3)-0.5*VoxSize));
elseif ZUp<=ZSide(2)&&YUp<=YSide(2)%
     Ytmp=ky*tUp+YSide(2);
     Ztmp=kz*tUp+ZSide(2);
     centerPoint(2)=point(2)+0.5*VoxSize-0.5*((point(2)+0.5*VoxSize)-Ytmp);
     centerPoint(3)=point(3)+0.5*VoxSize-0.5*((point(3)+0.5*VoxSize)-Ztmp);
     
 elseif ZUp>=ZSide(1)&&YUp<=YSide(2)%
     Ytmp=ky*tUp+YSide(2);
     %Ztmp=kz*tUp+ZSide(2);
     Ztmp=kz*tUp+ZSide(1);
     centerPoint(2)=point(2)+0.5*VoxSize-0.5*((point(2)+0.5*VoxSize)-Ytmp);
   %  centerPoint(3)=point(3)+0.5*VoxSize-0.5*((point(3)+0.5*VoxSize)-Ztmp);  
     centerPoint(3)=point(3)-0.5*VoxSize+0.5*(Ztmp-(point(3)-0.5*VoxSize));
     
     
     
  elseif YUp>=YSide(1)&&ZUp<=ZSide(2)%
     Ytmp=ky*tUp+YSide(1);
    centerPoint(2)=point(2)-0.5*VoxSize+0.5*(Ytmp-(point(2)-0.5*VoxSize));
        Ztmp=kz*tUp+ZSide(2);
    centerPoint(3)=point(3)+0.5*VoxSize-0.5*((point(3)+0.5*VoxSize)-Ztmp);
%    elseif YUp<=YSide(2)&&ZUp<=ZSide(2)%
%       Ytmp=ky*tUp+YSide(2);
%     centerPoint(2)=point(2)+0.5*VoxSize-0.5*((point(2)+0.5*VoxSize)-Ytmp);
%         Ztmp=kz*tUp+ZSide(1);
%     centerPoint(3)=point(3)-0.5*VoxSize+0.5*(Ztmp-(point(3)-0.5*VoxSize));
end
    
    
end

