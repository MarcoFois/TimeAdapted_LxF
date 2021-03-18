function F = flux(h,u,v)

g = 9.81;

F.hx = u;  F.hy = v;
F.ux = u.^2./(h+eps)+.5*g*h.^2; F.uy = u.*(v./(h+eps));
F.vx = v.*(u./(h+eps)); F.vy = v.^2./(h+eps)+.5*g*h.^2;

end
