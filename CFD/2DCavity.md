## *CFD Solution: 2D Lid-Driven Cavity*

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Use finite differences to approximate the Navier-Stokes equations, applied to a 2D lid-driven cavity.
                        The cavity has a horizontal velocity, u, of 1 at the top and vertical velocity, v = 0. The edges apply
                        a no-slip condition. Pressure, P, at the very bottom of the cavity is set equal to 0. The Reynolds 
                        number is 100.
                        <img align="center" width="100" height="100" src="/MATLAB/CFD/stagGrid.PNG">
                        
*Other problem-specific parameters:*
```MATLAB
Lx = 1;                         %size of x domain
Ly = 1;                         %size of y domain

nu = 0.01;                      %viscosity
rho = 1;                        %density

U_lid = 1;                      %lid velocity (top boundary)

dx = 0.05;                      %step size in x
dy = 0.05;                      %step size in y

C = 1;
dt = C*min([0.25*dx*dx/nu,dx/U_lid]); %time-step
```

![2DCavity](/MATLAB/CFD/2DCavity.png){:align = "center" height="50%" width="50%"}


Momentum and velocity are solved as the Burgers Equation using the conservative upwind scheme. The pressure poisson equation is 
solved using the SOR iterative method. P, u, and v are all solved using staggered grids.

![StagGrid](/MATLAB/CFD/stagGrid.PNG){:height="80%" width="80%"}

## Solution:

Below is a plot of the infinity norm of the residual returned from the SOR-pressure poisson solver. 
This shows great convergence of the pressure. Below that is plotted both U and V in the center of the
cavity compared to [trusted experimental data](https://www.sciencedirect.com/science/article/pii/0021999182900584).

![SORConverge](/MATLAB/CFD/SORConverge.png){:height="70%" width="70%"}

![UData](/MATLAB/CFD/UData.png){:height="50%" width="50%"}![VData](/MATLAB/CFD/VData.png){:height="50%" width="50%"}


*Example:*

The final results for u, v, and P are plotted individually below. As shown above, the resuts have converged to an accurate solution.

![UVelocity](/MATLAB/CFD/CavityU.png){:height="85%" width="85%"}

![VVelocity](/MATLAB/CFD/CavityV.png){:height="85%" width="85%"}

![Pressure](/MATLAB/CFD/CavityPressure.png){:height="85%" width="85%"}

## Implementation:
```MATLAB
%% PARAMETERS AND PROPERTIES

Lx = 1;                         %size of x domain
Ly = 1;                         %size of y domain

nu = 0.01;                      %viscosity
rho = 1;                        %density

U_lid = 1;                      %lid velocity (top boundary)

dx = 0.05;                      %step size in x
dy = 0.05;                      %step size in y

w = 1.6;                        %omega (SOR)

x_u = 0:dx:Lx;                  %discretized x for velocity in x direction
y_u = -dy/2:dy:Ly+dy/2;         %discretized y for velocity in x direction

x_v = -dx/2:dx:Lx+dx/2;         %discretized x for velocity in y direction
y_v = 0:dy:Ly;                  %discretized y for velocity in y direction

y_p = -dy/2:dy:Ly+dy/2;         %discretized y for pressure
x_p = -dx/2:dx:Lx+dx/2;         %discretized x for pressure

Tolfac = 1e-6;                  %Factor to determine PPE residual tolerance

gam = 0;                        %gamma for conservative donor cell approx.

C = 1;
dt = C*min([0.25*dx*dx/nu,dx/U_lid]); %time-step

t_f = 18;                       %end-time, after kinetic energy converged


imax = length(x_u);             %max i index inside physical domain
jmax = length(y_v);             %max j index inside physical domain

u = zeros(imax,jmax+1);         %initialize u and set boundaries
v = zeros(imax+1,jmax);         %initialize v and set boundaries

t=0:dt:t_f;                     %discretized time space
Nt = length(t);

for n = 1:Nt
    %boundary conditions for ghost nodes
    u(:,jmax+1) = 2*U_lid - u(:,jmax);
    u(:,1) = -u(:,2);
    
    v(1,:) = -v(2,:);
    v(imax+1,:) = -v(imax,:);
    
%% VELOCITY PREDICTOR (NON-LINEAR CONVECTION-DIFFUSION)
    
    F = u;                      %u-velocity predictor
    G = v;                      %v-velocity predictor
    
    for i = 2:imax - 1
        for j = 2:jmax
        %derivative terms in u-velocity
        d2udx2_p_d2udy2 = (1/(dx^2))*(u(i+1,j) - 2*u(i,j) + u(i-1,j))...
            + (1/(dy^2))*(u(i,j+1) - 2*u(i,j) + u(i,j-1));
        
        du2dx = (1/(4*dx))*((u(i,j) + u(i+1,j))^2-(u(i-1,j) + u(i,j))^2)...
            +(gam/(4*dx))*(abs(u(i,j)+u(i+1,j))*(u(i,j)-u(i+1,j))...
            -abs(u(i-1,j)+u(i,j))*(u(i-1,j)-u(i,j)));
        
        duvdy = (1/(4*dy))*((v(i,j)+v(i+1,j))*(u(i,j+1)+u(i,j))...
            -(v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j))) + (gam/(4*dy))...
            *(abs(v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1)) - abs(v(i,j-1)...
            +v(i+1,j-1))*(u(i,j-1)-u(i,j)));
        
        %u-velocity predictor
        F(i,j) = u(i,j) + dt*(nu*d2udx2_p_d2udy2 - du2dx - duvdy);    
        end
    end
    
    for i = 2:imax
        for j = 2:jmax - 1
        %derivative terms in v-velocity
        d2vdx2_p_d2vdy2 = (1/(dx^2))*(v(i+1,j) - 2*v(i,j) + v(i-1,j))...
            + (1/(dy^2))*(v(i,j+1) - 2*v(i,j) + v(i,j-1));
        
        dv2dy = (1/(4*dy))*((v(i,j) + v(i,j+1))^2-(v(i,j-1) + v(i,j))^2)...
            +(gam/(4*dy))*(abs(v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1))...
            -abs(v(i,j-1)+v(i,j))*(v(i,j-1)-v(i,j)));
        
        duvdx = (1/(4*dx))* ((v(i,j)+v(i+1,j))*(u(i,j+1)+u(i,j))...
            - (v(i-1,j)+v(i,j))*(u(i-1,j+1)+u(i-1,j))) + (gam/(4*dx))...
            *(abs(u(i,j+1)+u(i,j))*(v(i,j)-v(i+1,j)) - abs(u(i-1,j+1)...
            +u(i-1,j))*(v(i-1,j)-v(i,j)));
         
        %v-velocity predictor
        G(i,j) = v(i,j) + dt*(nu*d2vdx2_p_d2vdy2 - duvdx - dv2dy);
        end
    end 
    
%% PRESSURE POISSON

    %initialize pressure and source term
    S = zeros(imax+1,jmax+1);
    p = zeros(imax+1,jmax+1);
    
    %source term for PPE
    for i = 2:imax
        for j = 2:jmax
            
            S(i,j) = (rho/dt) * (((F(i,j) - F(i-1,j)) / (dx))...
                + ((G(i,j) - G(i,j-1)) / (dy)));
        end
    end
    
    %tolerance for iteration
    Tol = inf_norm(S)*Tolfac;
    
    %initialize residual and its norm to ensure iteration runs first time
    Res = p;
    Res_norm = 1e6;
    p_new = p;
    iter = 0;
    
    while (Res_norm >= Tol) && (iter < 1000)
        iter = iter + 1;
        for i = 2:imax
            for j = 2:jmax
                %boundary condition indicator functions
                eps_e = (i<imax);
                eps_w = (i>2);
                eps_n = (j<jmax);
                
                %SOR
                p_new(i,j) = p(i,j)*(1-w) + (w/(eps_w+eps_e+eps_n+1))...
                    *((eps_e*p(i+1,j)+eps_w*p_new(i-1,j))...
                    +(eps_n*p(i,j+1)+p_new(i,j-1))-dx*dx*S(i,j)); 
            end
        end
        
        %p is updated for next iteration
        p = p_new;
        
        for i = 2:imax
            for j = 2:jmax
                %boundary condition indicator functions
                eps_e = (i<imax);
                eps_w = (i>2);
                eps_n = (j<jmax);
                
                %Residual of calculated pressure
                Res(i,j) = ((eps_e*(p(i+1,j)-p(i,j))+eps_w*(p(i-1,j)-p(i,j))...
                    +eps_n*(p(i,j+1)-p(i,j))+(p(i,j-1)-p(i,j)))/(dx*dx))...
                    - S(i,j);
            end
        end
        
        Res_norm = inf_norm(Res);       %infinity norm of residual
        %store norm for later plotting
        plot_norm(iter) = Res_norm;
    end
    
%% SOLVE VELOCITY
    
    for i = 2:imax - 1
        for j = 2:jmax
            %u-velocity
            u(i,j) = F(i,j) - (dt/dx)*(p(i+1,j)-p(i,j));    
        end
    end
    
    for i = 2:imax
        for j = 2:jmax - 1
            %v-velocity
            v(i,j) = G(i,j) - (dt/dy)*(p(i,j+1)-p(i,j));
        end
    end

%% KINETIC ENERGY
    
    KE(n) = 0.5*rho*(sum(sum(u.^2)) + sum(sum(v.^2)));
end

%% PLOT RESULTS

%Kinetic Energy
figure
plot(t,KE)
xlabel("time (s)")
ylabel("KE")
title("Kinetic Energy")

%U Velocity
[X_U,Y_U] = meshgrid(x_u,y_u);
X_U=X_U';
Y_U=Y_U';
figure
surf(X_U,Y_U,u)
xlabel("x")
ylabel("y")
title("U")

%V Velocity
[X_V,Y_V] = meshgrid(x_v,y_v);
X_V=X_V';
Y_V=Y_V';
figure
surf(X_V,Y_V,v)
xlabel("x")
ylabel("y")
title("V")

%Pressure
[X_P,Y_P] = meshgrid(x_p,y_p);
X_P=X_P';
Y_P=Y_P';
figure
surf(X_P,Y_P,p)
xlabel("x")
ylabel("y")
title("Pressure")

%Compare to data
%u at x = 0.5
indx05 = find(x_u==0.5);
y_data=[1.0000 0.9766  0.9688 0.9609  0.9531 0.8516  0.7344 0.6172 0.5000...
    0.4531 0.2813 0.1719  0.1016 0.0703 0.0625 0.0547 0.0000];
u_data=[1.0000 0.8412 0.7887 0.7372 0.68717 0.2315 0.0033  -0.1364...
    -0.2058  -0.2109  -0.1566  -0.1015  -0.0643  -0.04775  -0.0419...
    -0.0371 0.0000];
figure
plot(y_data,u_data,'-o')
hold on
plot(y_u,u(indx05,:))
hold off
xlabel("y")
ylabel("u")
title("U: Results Compared to Data")
legend("data","results")

%v at y = 0.5
indy05 = find(y_v==0.5);
x_data=[1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5000...
    0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0.0000];
v_data=[0.0000 -0.05906  -0.0739 -0.0886 -0.10313 -0.16914 -0.22445...
    -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.1009...
    0.0923 0.0000];
figure
plot(x_data,v_data,'-o')
hold on
plot(x_v,v(:,indy05))
hold off
xlabel("x")
ylabel("v")
title("V: Results Compared to Data")
legend("data","results")

%Convergence of Pressure SOR
figure
semilogy(plot_norm)
xlabel("iterations")
ylabel("L_{inf}(Res)")
title("Convergence of Pressure SOR")

%% USER FUNCTIONS

%Error inf-norm
function f = inf_norm(A)
    f = max(max(abs(A)));
end
```
