function vnew = vstep (t, dt, x, dx, y, dy, h, u, v)

F = flux(h,u,v);

Flf_x = p0val (F.vx);
Flf_y = p0val (F.vy);
S   = Source(h, u, v,x,y);

g = 9.81;
% vel_star_x = max(norm(reshape(u+sqrt(g*h),1,numel(h)),inf), norm(reshape(u-sqrt(g*h),1,numel(h)),inf));
% vel_star_y = max(norm(reshape(v+sqrt(g*h),1,numel(h)),inf), norm(reshape(v-sqrt(g*h),1,numel(h)),inf));



vel_star_x = max(abs(u./(h+eps)+sqrt(g*h)), abs(u./(h+eps)-sqrt(g*h)));%max(abs(reshape(u+sqrt(g*h),1,numel(h))), abs(reshape(u-sqrt(g*h),1,numel(h))));
vel_star_y = max(abs(v./(h+eps)+sqrt(g*h)), abs(v./(h+eps)-sqrt(g*h)));%max(abs(reshape(v+sqrt(g*h),1,numel(h))), abs(reshape(v-sqrt(g*h),1,numel(h))));
vnew = astep (t, dt, x, dx, y, dy, v, Flf_x, Flf_y, S.v, vel_star_x, vel_star_y);
vnew(1,:)   = 0;
vnew(end,:) = 0;
end