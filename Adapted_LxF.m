clear all; close all; clc

% This script solves the "Sloshing problem" and the "dam break problem" in
% the context of shallow water equations using a time-adapted Lax-Friedrichs scheme.
% Copyright (c) 2021 - Marco Fois

%{ 
References

[1] M. Quecedo, M. Pastor. "A reappraisal of Taylor Galerkin algorithm for
drying-wetting areas in shallow water computations", Intl. Journal for Num.
Methods in Fluids 2002, 38:515-531.

[2] G.M. Porta, S. Perotto, F. Ballio. "A space-time adaptation scheme for
unsteady shallow water problems", Mathematics and Computers in simulation
2012, 82:2929-2950.

%}


% SPATIAL INFORMATIONS AND PHENOMENON -------------------------------------

fprintf('Which phenomenon do you want to simulate? \n')
fprintf('(1) Sloshing phenomenon \n');
fprintf('(2) Dam break problem \n')
n = input('');

switch n
    case 1
        fprintf('Sloshing problem \n');
        x = linspace(-400,400,50)';   dx = diff(x);
        y = linspace(0,80,20)';       dy = diff(y);
        fprintf('o----------------------------------------------o \n');
        fprintf('|                                              | \n');
        fprintf('|                                              | \n');
        fprintf('o----------------------------------------------o \n');
        fprintf('Rectangular domain: 800 x 80 m \n');
        
        [X,Y] = meshgrid(x,y);        [dX,dY] = meshgrid(dx,dy);

        Nele = (size(x,1)-1)*(size(y,1)-1);
        fprintf('Number of dofs: %g vertices. \n', numel(X));
        fprintf('Number of elements: %g squares. \n', Nele);
    case 2
        fprintf('Dam break problem \n');
        x = linspace(-5,5,33)';   dx = diff(x);
        y = linspace(-5,5,33)';   dy = diff(y);
        
        fprintf('o-------o \n');
        fprintf('|       | \n');
        fprintf('|       | \n');
        fprintf('o-------o \n');
        fprintf('Square domain: 5 x 5 m \n');
        
        [X,Y] = meshgrid(x,y);        [dX,dY] = meshgrid(dx,dy);

        Nele = (size(x,1)-1)*(size(y,1)-1);
        fprintf('Number of dofs: %g vertices. \n', numel(X));
        fprintf('Number of elements: %g squares. \n', Nele);
        
    otherwise
        error('Fool! You can choose only 1 or 2!');
              
end

% TEMPORAL INFORMATIONS ---------------------------------------------------

Time = 1.5;                       % Final time (s)
t = 0;                          % Initial time (s)
tol = 0.005 ;                     % Default tol for the adaptation

 
dt_min = 0.003;                  % Min dt allowed 
dt_max = 0.04;                     % Max dt allowed 
dt_k(1) = dt_min;               % First two time step
dt_k(2) = dt_min;
dtlist = [dt_k(1) dt_k(2)];
T = 0;


% INITIAL SOLUTIONS -------------------------------------------------------

switch n
    case 1
         h0 = 8 + sin(pi*X/800);
         u0 = zeros(size(X));
         v0 = zeros(size(X));

        h(:,:,1) = h0;        % Allocate the initial solution
        u(:,:,1) = u0;
        v(:,:,1) = v0;
        
        fprintf('The simulation ends after t = %g s. \n', Time);
        fprintf('Initial conditions: \n');
        fprintf('h(x,y) = 8 + sin(pi*x/800) \n');
        fprintf('u(x,y) = 0 \n');
        fprintf('v(x,y) = 0 \n');
        
    case 2
         h0 = ones(size(X)); %h0(X>=0 & X<=30) = 0.2;
         %h0(X>=-0.5 & X<=0.5 & Y>=-0.5 & Y<=0.5) = 2;
         h0(sqrt(X.^2+Y.^2)>=0 & sqrt(X.^2+Y.^2)<=0.5 ) = 2;
         u0 = zeros(size(X));
         v0 = zeros(size(X));

        h(:,:,1) = h0;        % Allocate the initial solution
        u(:,:,1) = u0;
        v(:,:,1) = v0;
        
        fprintf('The simulation ends after t = %g s. \n', Time);
        fprintf('Initial conditions: \n');
        fprintf('         { 2   for 0<= x^2 + y^2 <= 0.25 \n');
        fprintf('h(x,y) = {\n');
        fprintf('         { 1   elsewhere \n');
        fprintf('u(x,y) = 0 \n');
        fprintf('v(x,y) = 0 \n');

end

% FIRST THREE SOLUTIONS ---------------------------------------------------

for k = 1:2

    hnew(:,:) = hstep (t, dt_k(k), X, dX, Y, dY, h(:,:,k), u(:,:,k), v(:,:,k));
    unew(:,:) = ustep (t,dt_k(k), X, dX, Y, dY, h(:,:,k), u(:,:,k), v(:,:,k));
    vnew(:,:) = vstep (t, dt_k(k), X, dX, Y, dY, h(:,:,k), u(:,:,k), v(:,:,k));
    
    h(:,:,k+1) = hnew;
    u(:,:,k+1) = unew;
    v(:,:,k+1) = vnew;
    
    t = t + dt_k(k);
    T = [T t];
 
end

k = 3;
tic;

% ADAPTED SIMULATION ------------------------------------------------------

while t < Time
    
 % INTERPOLATIONS ---------------------------------------------------------
 
  % Time-derivative of the p1-interpolated solution between two time steps
  
        dh_t(:,:,k-1) = (h(:,:,k)-h(:,:,k-1))./(T(k)-T(k-1));
       
  % Time-derivative of the p2-interpolated solution between three time steps
 
        
        % Coefficients of the parabola h*(x,y,t)
        h1 = h(:,:,k-2)./((T(k-2)-T(k-1)).*(T(k-2)-T(k)));
        h2 = h(:,:,k-1)./((T(k-1)-T(k-2)).*(T(k-1)-T(k)));
        h3 = h(:,:,k)./((T(k)-T(k-2)).*(T(k)-T(k-1)));
        a(:,:) = h1+h2+h3; 
        b(:,:) = - (h1.*(T(k)+T(k-1)) + h2.*(T(k)+T(k-2)) + h3.*(T(k-1)+T(k-2)));
        
        
 % ERROR ESTIMATE nu_z* ---------------------------------------------------
      
        
  Nu_h(:,:,k-1) = (T(k)-T(k-1)).^2.*(4/3.*a.^2.*(T(k).^2+T(k).*T(k-1)+T(k-1).^2)+...
                     2.*a.*(b-dh_t(:,:,k-1)).*(T(k)+T(k-1))+(b-dh_t(:,:,k-1)).^2);
                

 
        
  nu_hmean(:,:,k-1) = p0val(Nu_h(:,:,k-1)); % Mean in each element 
  nu_htot(:,k-1) = sqrt(sum(nu_hmean(:,:,k-1),'all')); % Total estimator 
  rho_h(:,k-1) = (1/dt_k(k-1)).*sqrt(nu_htot(:,k-1)); % Computation of rho
  dt_k(k) = (tol./rho_h(:,k-1)).*sqrt(1/dt_k(k-1)); % Predict new dt

  
  
  if dt_k(k) > dt_max
      dt_k(k) = dt_max;
  elseif dt_k(k) < dt_min
      dt_k(k) = dt_min;
  end
  
        
  t = t + dt_k(k);
  dtlist = [dtlist dt_k(k)];
  T = [T t];
      
  
  % NEW SOLUTION ----------------------------------------------------------
  
    hnew = hstep (t, dt_k(k), X, dX, Y, dY, h(:,:,k), u(:,:,k), v(:,:,k));
    unew = ustep (t, dt_k(k), X, dX, Y, dY,  h(:,:,k), u(:,:,k), v(:,:,k));
    vnew = vstep (t, dt_k(k), X, dX, Y, dY,  h(:,:,k), u(:,:,k), v(:,:,k));
  
    k=k+1;
    
    h(:,:,k) = hnew;
    u(:,:,k) = unew;
    v(:,:,k) = vnew;
    

fprintf('Simulation time = %g, dt = %g \n', t, dt_k(k-1));

end
elapsed = toc;


fprintf('dt mean = %g \n', sum(dtlist)/numel(dtlist));
fprintf(' nu_htot mean = %g \n', sum(nu_htot)/numel(nu_htot));
fprintf('Elapsed time = %g  (s)\n', elapsed ); 

switch n
    
    case 1

 is = [1:50:numel(T)];
 for n = 1:numel(is)
       
    figure(1)
    surf (X, Y, h(:,:,is(n)))
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title(['Height h in time t = ' num2str(T(is(n))) ' [s]']);
    axis([-400 400 0 80 6.5 9.5])
    view(0,12)
    drawnow  
  
 end
 pause(1)
 close all;
 
 figure(2)
 subplot(2,1,1)
 plot(T(is),dtlist(is),'r-')
 legend('\Delta t')
 title('\Delta t step vs H^1(\Omega) error ')
 grid on
 subplot(2,1,2)
 plot(T(is),nu_htot(is),'b-')
 legend('\mu^T_{k-1} ')
 grid on
 
    case 2
        
 is = [1:1:numel(T)];
 for n = 1:numel(is)
       
    figure(1)
    surf (X, Y, h(:,:,is(n)))
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title(['Height h in time t = ' num2str(T(is(n))) ' [s]']);
    colorbar
    axis([-5 5 -5 5 0.5 2.5])
    view(45,45)
    drawnow  
  
 end
 pause(1)
 
 
 figure(2)
 subplot(2,1,1)
 plot(T(2:end),[ dtlist],'r-')
 legend('\Delta t')
 title('\Delta t step vs H^1(\Omega) error ')
 grid on
 subplot(2,1,2)
 plot(T(3:end),[ nu_htot],'b-')
 legend('\mu_{k-1}')
 grid on 
%   subplot(3,1,3)
%  semilogy(T(is),Nu_h(is),'b-')
%  legend('H^1 error')
%  grid on 

end



