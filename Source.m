function S = Source(h,u,v,x,y)

g = 9.81;
tau = .009;
%tau = 0;


% dZx = diff(Z,1,2)./diff(x,1,2);  dZy = diff(Z,1,1)./diff(y,1,1);
% dZx = [dZx(:,1),dZx];            dZy = [dZy,dZy(end,:)];   

S.h = 0*h;
S.u = -tau.*u;
S.v = -tau.*v;

end


