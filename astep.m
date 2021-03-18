function anew = astep (t, dt, x, dx, y, dy, a, Flf_x, Flf_y, S, vel_star_x, vel_star_y)
mass_x = .5 * ([dx, zeros(size(dx,1),1)] + [zeros(size(dx,1),1), dx]);
mass_x = [mass_x;mass_x(end,:)];
mass_y = .5 * ([dy; zeros(1,size(dy,2))] + [zeros(1,size(dy,2)); dy]);
mass_y = [mass_y,mass_y(:,end)];

mass = mass_x.*mass_y;



diffx = diff(vel_star_x.*a,1,2);
da_x = .5 * ( diffx(2:end,:) + diffx(1:end-1,:) );
diffy = diff(vel_star_y.*a,1,1);
da_y = .5 * ( diffy(:,2:end) + diffy(:,1:end-1) );

% viscosita' artificiale (ovvero flusso correttivo) fissata
% da = V .* diff(a);
   
% viscosita' artificiale (ovvero flusso correttivo) alla
% Lax-Friedrichs
% da_x = vel_star_x .* diffx;
% da_y = vel_star_y .* diffy;

% Fstar_x = Flf_x; 
% Fstar_y = Flf_y;

Fstar_x = Flf_x - da_x; 
Fstar_y = Flf_y - da_y;
 
Fstar_x = Fstar_x.*(dy/2);
Fstar_y = Fstar_y.*(dx/2); 
 
Fstar_x_ = zeros(size(a,1)+1, size(a,2)+1);
Fstar_x_(2:end-1,2:end-1) = Fstar_x;
Fstar_x_ = Fstar_x_(2:end,:)+Fstar_x_(1:end-1,:);
Fstar_x_ = -Fstar_x_(:,2:end) + Fstar_x_(:,1:end-1);
  
Fstar_y_ = zeros(size(a,1)+1, size(a,2)+1);
Fstar_y_(2:end-1,2:end-1) = Fstar_y;
Fstar_y_ = Fstar_y_(:,2:end)+Fstar_y_(:,1:end-1);
Fstar_y_ = -Fstar_y_(2:end,:) + Fstar_y_(1:end-1,:);
 
anew = a + dt * ( (Fstar_x_+Fstar_y_)./ mass + S );


end