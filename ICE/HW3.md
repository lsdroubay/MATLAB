**Analysis using Classical Lamination Theory**

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Calculate and plot the BMEP, Brake Power,Indicated Efficiencies, and Torque output of engine.

**Example:**

This example was done for a 4-stroke, 8-cylinder engine with bore and stroke of 10 cm and an assumed gross indicated efficiency (GIE) of 0.4.
The fuel used in this example is iso-octane, with its respective lower heating value and stoichiometric air-to-fuel ratio.

Some measured properties of the intake and exhaust are also included:
```MATLAB
P_exh = 1.1e5;  % measured exhaust pressure [Pa]

P_im = 0.58e5; % intake manifold pressure [Pa]
T_int = 330; % intake manifold temperature [K]
MW_int = 28.97; % molecular weight of intake charge [kg/kmol]
```

With these input parameters, the foloowing plots are generated with the MATLAB code following.
![BMEP](/MATLAB/ICE/BMEP.PNG)
![BMEP](https://user-images.githubusercontent.com/46492207/80048767-b7382800-84d6-11ea-8f79-55ee7b5366fc.png)


```MATLAB
close all;
clear all;
clc;

% Knowns
n_cyl = 8; % Number of cylinders
N = 2; % revolutions per cycle (4 -stroke --> n = 2, 2stroke --> n = 1)
bore = 100/1000; % cylinder bore [m]
stroke = 100/1000; % stroke [m]

AFR = 14.6; % air fuel ratio
LHV = 43e6; % fuels lower heating value [J/kg]

P_exh = 1.1e5;  % measured exhaust pressure [Pa]

P_im = 0.58e5; % intake manifold pressure [Pa]
T_int = 330; % intake manifold temperature [K]
MW_int = 28.97; % molecular weight of intake charge [kg/kmol]

GIE = 0.4; % gross indicated efficiency [-]

num_El = 1000;
RPM = linspace(1000,7000,num_El); % engine speed [rpm]

R_bar = 8314; % ideal gas constant [J/kmol-K]

% Correlations for FMEP and volumetric efficiency to engine speed
FMEP = (57+0.015*RPM+5.5e-6*RPM.^2) * 1000; % [Pa]
eta_vol = 2e-12*RPM.^3-5e-8*RPM.^2+0.0003*RPM+0.45;

%% SOLUTION
Vd = pi/4*bore^2*stroke*n_cyl; % engine displacement
rho_im = P_im/(R_bar/MW_int*T_int); % density of the intake air

%effective pressures, in units [Pa]
PMEP = P_im - P_exh;                  % pumping mean effectve pressure
IMEP_g = rho_im*eta_vol*GIE*LHV/AFR; % gross indicated mean effective pressure
IMEP_n = IMEP_g + PMEP;             % net indicated mean effective pressure
BMEP = IMEP_n - FMEP;                 % brake mean effective pressure

%efficiencies
eta_OCE = IMEP_n./IMEP_g;           % open cycle efficiency
eta_ME = BMEP./IMEP_n;              % mechanical efficiency

% brake engine parameters
W_b = BMEP*Vd;                      % brake work
W_dot_brake = W_b.*RPM/(N*60);      % brake power
T_b = W_dot_brake*60./(2*pi*RPM);   % brake torque
BTE = GIE.*eta_ME.*eta_OCE;         % brake thermal efficiency

%% PLOTTING
% BMEP
figure;
plot(RPM,BMEP/1e5,'k','linewidth',2)
xlabel('Engine Speed [rev/min]')
ylabel('BMEP [bar]')

% Brake Power
figure;
plot(RPM,W_dot_brake/1000,'k','linewidth',2)
xlabel('Engine Speed [rpm]')
ylabel('Brake Power [kW]')

%Efficiencies
figure
plot(RPM,BTE*100,'k','linewidth',2)
hold all;
plot(RPM,BTE*100./eta_ME,'r','linewidth',2)
plot(RPM,BTE*100./(eta_ME.*eta_OCE),'b','linewidth',2)
legend({'BTE';'NIE';'GIE'},'Location','Best')
xlabel('Engine Speed [rpm]')
ylabel('Efficiency [%]')

% Torque
figure;
plot(RPM,T_b,'k','linewidth',2)
xlabel('Engine Speed [rpm]')
ylabel('Torque [N-m]')
```
