function [ur,vr,thetad]=my_rotate_pres(u,v,p,thetad)
%function [ur,[vr],thetad]=rotate(u,[v])
%           Rotates into direction of maximum variance if no theta
%           specified, if theta srotate sin to that direction

%     u is vector velocity, ur is rotated vector velocity
%      direction of theta is counterclockwise, based on E vel 1st, N vel second
if nargin==1;
    v=imagg(u);
    u=real(u);
        maxvar=1;

end
if nargin==2|3;
            maxvar=1;

end
if nargin==4;
    maxvar=0;
end

if any( isnan(u(:)))|(isnan(v(:)) )
    disp('NaNs found')
    ur=ones(size(u))*nan;
    vr=ur;
    g=find( (~isnan(u(:)))&(~isnan(v(:))) );
    if length(g)>1;
        u=u(g);
        v=v(g);
    end
else
    ur=u;
    vr=u;
    g=1:length(u(:));
end
if maxvar
    if length(g)>1
        up=cov([u(:) p(:)]);
                vp=cov([v(:) p(:)]);

        %cv=cov([u(:) v(:)]);
       theta=atan2(vp(1,2),up(1,2));
        
        thetad=theta*180./pi;

        
    else;
        thetad=nan;
    end
end

vr(g)=-u*sind(thetad)+v*cosd(thetad);
ur(g)=u*cosd(thetad)+v*sind(thetad);

if nargin==1
    ur=ur+sqrt(-1)*vr;
    vr=thetad;
    thetad=[];
end

