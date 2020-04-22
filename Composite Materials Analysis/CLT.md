**Analysis using Classical Lamination Theory**

**Author:** Landon Droubay

**Language:** MATLAB


**Description/Purpose:** Returns the needed properties and values for composite material laminate. Uses created [LaminaClass](/LaminaClass) class.

**Output:** This displays the material properties, laminate properties, Qb matrix, loads and moments, 
            ABD matrix, strains and curvatures, normal and shear stresses, normal and shear strain, 
            shear and strain in   , average/smear properties, and the failure factors/criteria.

**Usage/Example:**

For a symmetric carbon fiber 5-layer laminate of angles: 30,-30,0,-30,30. Layer thickness is 1.5e-4 m

```
LaminateAngles = [30,-30,0,-30,30];
symmetric = true;


N=numel(LaminateAngles);            %Number of layers

for i=1:N
    k(i) = LaminaClass;
end


%initialize angles, materials, and thicknesses for each layer

for i = 1:N
    k(i).theta = LaminateAngles(i);
end

for i = 1:N
    k(i).material = 1;          %1 = graphite, 2 = glass, 3 = other
end

for i = 1:N
    k(i).thickness = 1.5*(10^-4);       %in meters
end

for i = 1:N
    k(i) = props(k(i));                 %Function Below
end

%Display Material Properties
DispProp(k);                            %Function Below

%z values for each layer
z = zeros(1,N+1);

EvenOdd = mod(N,2);
if EvenOdd == 0    
    for i = (N/2):-1:1
        z(i) = z(i+1) - k(i).thickness;
    end

    for i = ((N/2)+2):(N+1)
        z(i) = z(i-1) + k(i-1).thickness;
    end
else
    mid =(N-1)/2 ;
    i = mid;
    z(i+1) = -k(i).thickness/2;
    z(i+2) = k(i+1).thickness/2;
    
    for i = mid:-1:1
        z(i) = z(i+1) - k(i).thickness;
    end
    
    for i = (mid + 2): (N+1)
        z(i) = z(i-1) + k(i-1).thickness;
    end
end




%Array of Qbar matrices for each layer

for i = 1:N
    Qb(:,:,i) = Qbar(k(i));             %Function Below
    
    %Print Results for Each Layer
    fprintf('Lamina %i\ttheta = %i\tt = %f\tmaterial = %i\n',i,k(i).theta,k(i).thickness,k(i).material)
    fprintf('Qb(%i)= \t\t\t\t\tin Pa\n',i)
    disp(Qb(:,:,i))
end


%Loads and Moments
NM = [10000;-60000;0;1;-3;0];

%Thermal Loads
NMTherm = ThermLoads(Qb,k,z);

deltaT = -120;
NMT = NMTherm*deltaT;

%Print loads
fprintf('Nx=\t%d\tNy=\t%d\tNxy=\t%d\nMx=\t%d\tMy=\t%d\tMxy=\t%d\n',...
    NM(1),NM(2),NM(3),NM(4),NM(5),NM(6))
fprintf('NxT=\t%d\tNyT=\t%d\tNxyT=\t%d\nMxT=\t%d\tMyT=\t%d\tMxyT=\t%d\ndeltaT=%f\n',...
    NMT(1),NMT(2),NMT(3),NMT(4),NMT(5),NMT(6),deltaT)

%Combine Loads
NM = NM + NMT;

%Given midplane strains and curvatures
EK0 = [0;0;0;0;0;0];

%ABD Matrix for Laminate
ABDmat = ABD(Qb,z);                     %Function Below

%abd Matrix
abd = inv(ABDmat);
disp('abd = ')
disp(abd)



EK0e = abd*NM;         %Strain calculated using known loads and abd matrix         
NMe = ABDmat*EK0;      %Loads calculated using known strains and ABD matrix

%Display Found Loads
fprintf('Strains & Curvatures\ne0x=\t%d\te0y=\t%d\tg0xy=\t%d\nK0x=\t%d\tK0y=\t%d\tK0xy=\t%d\n',...
    EK0e(1),EK0e(2),EK0e(3),EK0e(4),EK0e(5),EK0e(6))

%Calculate and display stresses and strains; off-axis and on-axis
sig12 = getStressStrain(Qb,EK0e,k,z,deltaT);

%Display Smear Properties if Symmetric
if symmetric == true
    SmearProps(z,abd,NMTherm,k,sig12,deltaT,EK0e)
end

%Check and Display Failure Criterion
Failure(Qb,EK0e,k,z);




%FUNCTIONS
%__________________________________________________________________________
function p = props(k)

mat = k.material;
    %Case for different materials
    switch mat
        case 1                          %Graphite
            E(1) = 155;
            E(2) = 12.10;
            E(3) = E(2);

            V(2,3) = 0.458;
            V(1,3) = 0.248;
            V(1,2) = 0.248;

            G(1) = 3.2;  %G23
            G(2) = 4.4;  %G13
            G(3) = 4.4;  %G12

            a1 = -0.01800 * (10^-6);
            a2 = 24.3 * (10^-6);
            a3 = a2;
            
            sigmaT1 = 1500 * (10^6);
            sigmaC1 = -1250 * (10^6);
            sigmaT2 = 50 * (10^6);
            sigmaC2 = -200 * (10^6);
            tauF12 = 100 * (10^6);

        case 2                          %Glass
            E(1) = 50;
            E(2) = 15.20;
            E(3) = E(2);

            V(2,3) = 0.428;
            V(1,3) = 0.254;
            V(1,2) = 0.254;

            G(1) = 3.28;  %G23
            G(2) = 4.7;  %G13
            G(3) = 4.7;  %G12

            a1 = 6.34 * (10^-6);
            a2 = 23.3 * (10^-6);
            a3 = a2;
            
            sigmaT1 = 1000 * (10^6);
            sigmaC1 = -6000 * (10^6);
            sigmaT2 = 30 * (10^6);
            sigmaC2 = -120 * (10^6);
            tauF12 = 70 * (10^6);

        case 3                          %Other
            
            %Prompt to Get Values
            
            %get Es
            prompt = 'What is E1? ';
            E(1) = input(prompt);

            prompt = 'What is E2? ';
            E(2) = input(prompt);

            prompt = 'What is E3? ';
            E(3) = input(prompt);

            %get Poissons
            V = zeros(3,3);

            for i = 1:3
                for j = 1:3
                    if i~=j
                        prompt = sprintf('Do you know v%d%d? Y/N ', i, j);
                        str = input(prompt,'s');
            if str == 'Y'
                 prompt = sprintf('What is v%d%d? ',i,j);
                 V(i,j) = input(prompt);
            end
                    end
                end
            end


            %get Gs
            prompt = 'What is G1? ';
            G(1) = input(prompt);

            prompt = 'What is G2? ';
            G(2) = input(prompt);

            prompt = 'What is G3? ';
            G(3) = input(prompt);
        otherwise
            disp('other value')
    end

E = E * (10^9);
G = G * (10^9);


    %find unknown Poissons
    for i = 1:3
        for j = 1:3
            if i~=j
                if V(i,j) ~= 0

                    V(j,i) = (V(i,j) * E(j)) / E(i);

                end
            end
        end
    end

    %Set up Compliance Matrix
    S = zeros(6,6);

    for i = 1:3

        S(i,i) = 1 / E(i);

    end

    for i = 1:3
        for j = 1:3
            if i~=j

                S(i,j) = -(V(i,j)) / E(i);

            end
        end
    end

    for i = 4:6

        j = i-3;
       S(i,i) = 1 / G(j);

    end

    %Stiffness Matrix
    C = inv(S);

    %Change S to (TPa)^-1
    %S = S*1000;

    %Reduced Compliance Matrix
    for i =1:2
        for j=1:2
            Sreduced(i,j) = S(i,j);
        end
    end
    Sreduced(1,3) = S(1,6);
    Sreduced(3,1) = S(1,6);
    Sreduced(2,3) = S(2,6);
    Sreduced(3,2) = S(2,6);
    Sreduced(3,3) = S(6,6);

    %Reduced Stiffness Matrix
    Q = zeros(3,3);
    Q(1,1) = C(1,1) -((C(1,2)*C(2,1)) / C(2,2));
    Q(1,2) = C(1,2) -((C(1,2)*C(2,3)) / C(2,2));
    Q(2,1) = Q(1,2);
    Q(2,2) = C(2,2) -((C(2,3)*C(2,3)) / C(2,2));
    Q(3,3) = C(6,6);
    
    k.E1 = E(1);
    k.E2 = E(2);
    k.v12 = V(1,2);
    k.v21 = V(2,1);
    k.G12 = G(3);
    k.Q = Q;
    k.S = S;
    k.a(1,1) = a1;
    k.a(2,1) = a2;
    k.a(3,1) = a3;
    k.axy(1,1) = (cosd(k.theta)^2)*a1 + (sind(k.theta)^2)*a2;
    k.axy(2,1) = (sind(k.theta)^2)*a1 + (cosd(k.theta)^2)*a2;
    k.axy(3,1) = 2*cosd(k.theta)*sind(k.theta)*(a1-a2);
    k.sigmaT1 = sigmaT1;
    k.sigmaC1 = sigmaC1;
    k.sigmaT2 = sigmaT2;
    k.sigmaC2 = sigmaC2;
    k.tauF12 = tauF12;
    
    p = k;
end

function DispProp(kk)
    
   num = numel(kk);
   mats = zeros(1,num);
   
   for i = 1:num
       if(kk(i).material ~= mats(1:num))
          fprintf('Material = %i\nE1 = %e\tE2 = %e\tG12 = %e\n'...
              ,kk(i).material,kk(i).E1,kk(i).E2,kk(i).G12)
          fprintf('v12 = %e\tv21 = %e\n',kk(i).v12,kk(i).v21)
          fprintf('Q(%i) = \t\t\t\t\tin Pa\n',kk(i).material)
          disp(kk(i).Q)
          
           
       end
       
       mats(i) = kk(i).material;
       Lam(i) = i;
       Thick(i) = kk(i).thickness;
       Angle(i) = kk(i).theta;
       
       Material = mats.';
       Lamina = Lam.';
       Thickness = Thick.';
       Orientation = Angle.';
   end
   fprintf('\n\n')
   T = table(Lamina, Material, Thickness, Orientation);
   disp(T)
   fprintf('\n')
end


function QQ = Qbar(k)




    
    theta = k.theta;
    S = k.S;
    Q = k.Q;
    

    m = cosd(theta);            %m for transformations (calculated for every theta)
    n = sind(theta);            %n for transformations (calculated for every theta)
    
    
    SS = zeros(3,3);
    %S-bar matrix for all theta  
        SS(1,1) = (S(1,1) * (m^4)) + ((2*S(1,2) + S(6,6)) * (m^2) * (n^2)) + (S(2,2)*(n^4));

        SS(1,2) = ((S(1,1) + S(2,2) - S(6,6)) * (m^2) * (n^2)) + (S(1,2)...
        * ((n^4) + (m^4)));

        SS(2,1) = SS(1,2);

        SS(2,2) = (S(1,1) * (n^4)) + ((2*S(1,2) + S(6,6)) * (m^2) * (n^2)) + (S(2,2)*(m^4));

        SS(1,3) = (((2*S(1,1)) - (2*S(1,2)) - (S(6,6))) * n * m^3) + ...
        (((2*S(1,2)) - (2*S(2,2)) + (S(6,6))) * m * n^3);

        SS(3,1) = SS(1,3);

        SS(2,3) = (((2*S(1,1)) - (2*S(1,2)) - (S(6,6))) * m * n^3) + ...
        (((2*S(1,2)) - (2*S(2,2)) + (S(6,6))) * n * m^3);

        SS(3,2) = SS(2,3);

        SS(3,3) = (2*((2*S(1,1)) + (2 * S(2,2)) - (4*S(1,2)) - S(6,6)) ...
        * (m^2) * (n^2)) + (S(6,6) * ((n^4) + (m^4)));



    %Q-bar matrix for all theta

        Qb(1,1) = (Q(1,1) * (m^4)) + (((2*Q(1,2)) + (4*Q(3,3))) * (m^2) * (n^2)) + (Q(2,2)*(n^4));

        Qb(1,2) = ((Q(1,1) + Q(2,2) - (4*Q(3,3))) * (m^2) * (n^2)) + (Q(1,2)...
        * ((n^4) + (m^4)));

        Qb(2,1) = Qb(1,2);

        Qb(2,2) = (Q(1,1) * (n^4)) + ((2*Q(1,2) + (4*Q(3,3))) * (m^2) * (n^2)) + (Q(2,2)*(m^4));

        Qb(1,3) = ((Q(1,1) - Q(1,2) - (2*Q(3,3))) * (m^3) * n) + ((Q(1,2) - Q(2,2) + (2*Q(3,3))) * (n^3) * m);

        Qb(3,1) = Qb(1,3);

        Qb(2,3) = ((Q(1,1) - Q(1,2) - (2*Q(3,3))) * (n^3) * m) + ((Q(1,2) - Q(2,2) + (2*Q(3,3))) * (m^3) * n);

        Qb(3,2) = Qb(2,3);

        Qb(3,3) = (((Q(1,1)) + (Q(2,2)) - (2*Q(1,2)) - (2*Q(3,3))) ...
        * (m^2) * (n^2)) + (Q(3,3) * ((n^4) + (m^4)));


    QQ = Qb;
end


function ABDmat = ABD(Qb,z)

layers = (numel(z)) - 1;

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

%A Matrix
for k = 1:layers
    A = A + (Qb(:,:,k) * (z(k+1) - z(k)));
end


%B Matrix
for k = 1:layers
    B = B + (Qb(:,:,k) * (z(k+1)^2 - z(k)^2));
end
B = B/2;


%D Matrix
for k = 1:layers
    D = D + (Qb(:,:,k) * (z(k+1)^3 - z(k)^3));
end

D = D/3;


%Full ABD matrix
A_B_D = zeros(6,6);
for i=1:3
    for j =1:3
        A_B_D(i,j) = A(i,j);
    end
end

for i=1:3
    for j = 4:6
        A_B_D(i,j) = B(i,j-3);
    end
end

for i=4:6
    for j = 1:3
        A_B_D(i,j) = B(i-3,j);
    end
end

for i = 4:6
    for j = 4:6
    A_B_D(i,j) = D(i-3,j-3);
    end
end

%display A, B, and D
A
B
D

ABDmat = A_B_D;

end


function SmearProps(z,abd,NMTherm,k,sig12,deltaT,EK0e)

H = z(1) * -2;

%Due to Sigma Xbar
EXbar = 1 / (H * abd(1,1));
vXYbar = -abd(1,2) / abd(1,1);
NxyxBar = abd(1,3) / abd(1,1);

%Due to Sigma Ybar
EYbar = 1 / (H * abd(2,2));
vYXbar = -abd(1,2) / abd(2,2);
NxyyBar = abd(2,3) / abd(2,2);

%Due to Tau XYbar
GXYbar = 1 / (H * abd(3,3));
NxxyBar = -abd(1,3) / abd(3,3);
NyxyBar = abd(2,3) / abd(3,3);

%Through Thickness
ezbar = vtable(k,sig12);
vxzbar = -ezbar/EK0e(1);
vyzbar = -ezbar/EK0e(2);

%Thermal Properties
axbar = abd(1,1) * NMTherm(1) + abd(1,2)*NMTherm(2) + abd(1,3)*NMTherm(3);
aybar = abd(2,1) * NMTherm(1) + abd(2,2)*NMTherm(2) + abd(2,3)*NMTherm(3);
axybar = abd(3,1) * NMTherm(1) + abd(3,2)*NMTherm(2) + abd(3,3)*NMTherm(3);

ezbarT = atable(k,sig12,deltaT);
azbar = ezbarT / deltaT;

%Display Smear Properties
fprintf('\nSmear Properties\nEx=\t%d\tEy=\t%d\tGxy=\t%d\nvxy=\t%d\tvyx=\t%d\n',...
    EXbar,EYbar,GXYbar,vXYbar,vYXbar)
fprintf('vxz=\t%d\tvyz=\t%d\n',vxzbar,vyzbar)
fprintf('Nxyx=\t%d\tNxyy=\t%d\tNxxy=\t%d\tNyxy=\t%d\n'...
    ,NxyxBar,NxyyBar,NxxyBar,NyxyBar)
fprintf('ax=\t%d\tay=\t%d\taxy=\t%d\taz=\t%d\n',axbar,aybar,axybar,azbar)

end

function T = Trans(theta)

T = [cosd(theta)^2, sind(theta)^2,2*sind(theta)*cosd(theta);sind(theta)^2,cosd(theta)^2,-2*cosd(theta)*sind(theta);-cosd(theta)*sind(theta),cosd(theta)*sind(theta),cosd(theta)^2-sind(theta)^2];

end

function Tr = TransStrain(theta)

T = [cosd(theta)^2, sind(theta)^2,2*sind(theta)*cosd(theta);sind(theta)^2,cosd(theta)^2,-2*cosd(theta)*sind(theta);-cosd(theta)*sind(theta),cosd(theta)*sind(theta),cosd(theta)^2-sind(theta)^2];
R(1,1) = 1;
R(2,2) = 1;
R(3,3) = 2;

R1(1,1) = 1;
R1(2,2) = 1;
R1(3,3) = 0.5;

Tr = R*T*R1;
end

function therm = ThermLoads(Qb,k,z)
N = numel(k);

NTherm = [0;0;0];

for i = 1:N
    NTherm = NTherm + (Qb(:,:,i) * k(i).thickness * k(i).axy);
end

MTherm = [0;0;0];

for i = 1:N
    MTherm = MTherm + ((Qb(:,:,i) * k(i).axy) * 0.5 *(z(i+1)^2 - z(i)^2));
end

for i =1:3
    therm(i,1) = NTherm(i);
end

for i = 4:6
    therm(i,1) = MTherm(i-3);
end

end

function sig12 = getStressStrain(Qb,ek,k,z,T)
        N = numel(z) - 1;

        E = [ek(1);ek(2);ek(3)];
        K = [ek(4);ek(5);ek(6)];
        
        for i = 1:N%+1
            eps(:,:,i) = E + z(i) * K - (T*k(i).axy) ;
        end
eps(:,:,N+1) = E + K*z(N+1) - (T*k(N).axy);
        fprintf('\nLamina\tLocation\tsigmax\t\t\tsigmay\t\t\ttauxy\n')
        for i = 1:N
            pos = 'top';
            Nsig = Qb(:,:,i) * eps(:,:,i);
            Nsig12 = Trans(k(i).theta) * Nsig;
            fprintf('%i\t\t%s\t\t%d\t\t%d\t\t%d\n',i,pos,Nsig(1),Nsig(2),Nsig(3))
            pos = 'bottom';
            Nsig = Qb(:,:,i) * eps(:,:,i+1);
            fprintf('%i\t\t%s\t%d\t\t%d\t\t%d\n',i,pos,Nsig(1),Nsig(2),Nsig(3))
            
            for l = 1:3
    
                sig(l,1,i) = Nsig(l);
                sig12(l,1,i) = Nsig12(l);
            end
            
            %fprintf
        end
        
        fprintf('\nLamina\tLocation\tepsx\t\t\tepsy\t\t\tgamxy\n')
        for i=1:N
            pos = 'top';
            fprintf('%i\t\t%s\t\t%d\t\t%d\t\t%d\n',i,pos,eps(1,1,i),eps(2,1,i),eps(3,1,i))
            pos = 'bottom';
            fprintf('%i\t\t%s\t%d\t\t%d\t\t%d\n',i,pos,eps(1,1,i+1),eps(2,1,i+1),eps(3,1,i+1))
        end
        
         fprintf('\nLamina\tLocation\tsigma1\t\t\tsigma2\t\t\ttau12\n')
         for i = 1:N
            pos = 'top';
            Nsig = Qb(:,:,i) * eps(:,:,i);
            Nsig12 = Trans(k(i).theta) * Nsig;
            fprintf('%i\t\t%s\t\t%d\t\t%d\t\t%d\n',i,pos,Nsig12(1),Nsig12(2),Nsig12(3))
            pos = 'bottom';
            Nsig = Qb(:,:,i) * eps(:,:,i+1);
            Nsig12 = Trans(k(i).theta) * Nsig;
            fprintf('%i\t\t%s\t%d\t\t%d\t\t%d\n',i,pos,Nsig12(1),Nsig12(2),Nsig12(3))
         end
           
         fprintf('\nLamina\tLocation\teps1\t\t\teps2\t\t\tgam12\n')
         for i = 1:N
             eps = E + z(i) * K;
             eps12 = TransStrain(k(i).theta) * eps;
             pos = 'top';
             fprintf('%i\t\t%s\t\t%d\t\t%d\t\t%d\n',i,pos,eps12(1),eps12(2),eps12(3))
             pos = 'bottom';
             eps = E + z(i+1) * K;
             eps12 = TransStrain(k(i).theta) * eps;
             fprintf('%i\t\t%s\t%d\t\t%d\t\t%d\n',i,pos,eps12(1),eps12(2),eps12(3))
        
         end

end

function Failure(Qb,ek,k,z)

N = numel(z) - 1;

        E = [ek(1);ek(2);ek(3)];
        K = [ek(4);ek(5);ek(6)];
        
        
         fprintf('\nLamina\tLocation\tfac1\t\t\tfac2\t\t\t\tfac12\t\t\tTsai Wu\n')
         for i = 1:N
            pos = 'top';
            Nsig = Qb(:,:,i) * (E + z(i+1)*K);
            Nsig12 = Trans(k(i).theta) * Nsig;
            TW = TsaiWu(Nsig12,k(i));
            
            if Nsig12(1) >= 0 
               fac1 = Nsig12(1) /  k(i).sigmaT1;
            else
               fac1 = Nsig12(1) /  k(i).sigmaC1;
            end
            
            if Nsig12(2) >= 0 
               fac2 = Nsig12(2) /  k(i).sigmaT2;
            else
               fac2 = Nsig12(2) /  k(i).sigmaC2;
            end
            
            fac12 = Nsig12(3) / k(i).tauF12;
            
            fprintf('%i\t\t%s\t\t%d\t\t%d\t\t%d\t\t%d\n',i,pos,fac1,fac2,fac12,TW)
            pos = 'bottom';
            Nsig = Qb(:,:,i) * (E + z(i+1)*K);
            Nsig12 = Trans(k(i).theta) * Nsig;
            TW = TsaiWu(Nsig12,k(i));
            
            if Nsig12(1) >= 0 
               fac1 = Nsig12(1) /  k(i).sigmaT1;
            else
               fac1 = Nsig12(1) /  k(i).sigmaC1;
            end
            
            if Nsig12(2) >= 0 
               fac2 = Nsig12(2) /  k(i).sigmaT2;
            else
               fac2 = Nsig12(2) /  k(i).sigmaC2;
            end
            
            fac12 = Nsig12(3) / k(i).tauF12;
            
            fprintf('%i\t\t%s\t%d\t\t%d\t\t%d\t\t%d\n',i,pos,fac1,fac2,fac12,TW)
         end
           
      



end

function facTW = TsaiWu(sig,k)

F1 = (1/k.sigmaT1) + (1/k.sigmaC1);
F2 = (1/k.sigmaT2) + (1/k.sigmaC2);
F11 = -1/(k.sigmaT1*k.sigmaC1);
F22 = -1/(k.sigmaT2*k.sigmaC2);
F66 = (1/k.tauF12)^2;
F12 = -0.5 * sqrt(F11*F22);

facTW = F1*sig(1) + F2*sig(2) + F11*(sig(1)^2) + F22*(sig(2)^2) + F66*(sig(3)^2) + 2*F12*sig(1)*sig(2);

end

function ez = vtable(k,sig)
N= numel(k);
HH = N*k(1).thickness;
H = 0;
for i = 1:N
   e3(i) = k(i).S(1,3)*sig(1,1,i) + k(i).S(2,3)*sig(2,1,i);
   h(i) = e3(i)*k(i).thickness; 
   H = H + h(i);
end
ez = H/HH;
end

function ez = atable(k,sig,T)
N= numel(k);
HH = N*k(1).thickness;
H = 0;
for i = 1:N
   e3(i) = k(i).S(1,3)*sig(1,1,i) + k(i).S(2,3)*sig(2,1,i) + k(i).a(3)*T;
   h(i) = e3(i)*k(i).thickness; 
   H = H + h(i);
end

ez = H/HH;
end
```

The reults for this example:

```
Material = 1
E1 = 1.550000e+11	E2 = 1.210000e+10	G12 = 4.400000e+09
v12 = 2.480000e-01	v21 = 1.936000e-02
Q(1) = 					in Pa
   1.0e+11 *

    1.5575    0.0302         0
    0.0302    0.1216         0
         0         0    0.0440



    Lamina    Material    Thickness    Orientation
    ______    ________    _________    ___________

      1          1         0.00015          30    
      2          1         0.00015         -30    
      3          1         0.00015           0    
      4          1         0.00015         -30    
      5          1         0.00015          30    


Lamina 1	theta = 30	t = 0.000150	material = 1
Qb(1)= 					in Pa
   1.0e+10 *

    9.2799    3.0067    4.6706
    3.0067    2.1004    1.5470
    4.6706    1.5470    3.1452

Lamina 2	theta = -30	t = 0.000150	material = 1
Qb(2)= 					in Pa
   1.0e+10 *

    9.2799    3.0067   -4.6706
    3.0067    2.1004   -1.5470
   -4.6706   -1.5470    3.1452

Lamina 3	theta = 0	t = 0.000150	material = 1
Qb(3)= 					in Pa
   1.0e+11 *

    1.5575    0.0302         0
    0.0302    0.1216         0
         0         0    0.0440

Lamina 4	theta = -30	t = 0.000150	material = 1
Qb(4)= 					in Pa
   1.0e+10 *

    9.2799    3.0067   -4.6706
    3.0067    2.1004   -1.5470
   -4.6706   -1.5470    3.1452

Lamina 5	theta = 30	t = 0.000150	material = 1
Qb(5)= 					in Pa
   1.0e+10 *

    9.2799    3.0067    4.6706
    3.0067    2.1004    1.5470
    4.6706    1.5470    3.1452

Nx=	10000	Ny=	-60000	Nxy=	0
Mx=	1	My=	-3	Mxy=	0
NxT=	-1.039078e+04	NyT=	-2.253681e+04	NxyT=	0
MxT=	1.040834e-16	MyT=	0	MxyT=	0
deltaT=-120.000000

A =

   1.0e+07 *

    7.9041    1.8492   -0.0000
    1.8492    1.4426   -0.0000
   -0.0000   -0.0000    1.9531


B =

   1.0e-12 *

   -0.9095         0         0
         0   -0.1137         0
         0         0         0


D =

    3.2802    1.0494    0.9458
    1.0494    0.7359    0.3133
    0.9458    0.3133    1.0981

abd = 
    0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000
   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000
    0.0000    0.0000    0.0000    0.0000   -0.0000   -0.0000
    0.0000   -0.0000    0.0000    0.6555   -0.7904   -0.3391
   -0.0000    0.0000   -0.0000   -0.7904    2.4997   -0.0323
   -0.0000    0.0000   -0.0000   -0.3391   -0.0323    1.2119

Strains & Curvatures
e0x=	1.904899e-03	e0y=	-8.163146e-03	g0xy=	-2.591886e-20
K0x=	3.026784e+00	K0y=	-8.289496e+00	K0xy=	-2.421946e-01

Lamina	Location	sigmax		sigmay		tauxy
1	top	-6.108849e+07	-5.291506e+07	-5.106697e+07
1	bottom	1.780334e+08	1.224712e+07	1.087299e+08
2	top	-6.312941e+07	-6.762910e+07	5.366737e+07
2	bottom	1.558455e+07	-4.704976e+07	-6.141659e+06
3	top	2.470443e+08	-5.118497e+07	7.992422e+04
3	bottom	4.254317e+08	-7.160426e+07	1.103976e+07
4	top	-5.024323e+07	-9.143714e+07	4.744203e+07
4	bottom	1.922723e+08	-2.515095e+07	-1.146400e+08
5	top	-4.889054e+07	-1.050272e+08	-4.775720e+07
5	bottom	-4.584105e+07	-1.180552e+08	-4.692976e+07

Lamina	Location	epsx		epsy		gamxy
1	top	1.497235e-03	-2.868125e-03	-2.436378e-03
1	bottom	1.951252e-03	-4.111550e-03	2.581694e-03
2	top	1.951252e-03	-4.111550e-03	2.581694e-03
2	bottom	1.675730e-03	-4.625434e-03	1.816459e-05
3	top	1.675730e-03	-4.625434e-03	1.816459e-05
3	bottom	2.859288e-03	-6.598398e-03	2.509036e-03
4	top	2.859288e-03	-6.598398e-03	2.509036e-03
4	bottom	3.313305e-03	-7.841823e-03	-2.581694e-03
5	top	3.313305e-03	-7.841823e-03	-2.581694e-03
5	bottom	3.767323e-03	-9.085247e-03	-2.618024e-03

Lamina	Location	sigma1		sigma2		tau12
1	top	-1.032704e+08	-1.073313e+07	-2.199429e+07
1	bottom	2.307497e+08	-4.046913e+07	-1.742265e+07
2	top	-1.107316e+08	-2.002687e+07	2.878211e+07
2	bottom	5.244802e+06	-3.671002e+07	2.405062e+07
3	top	2.470443e+08	-5.118497e+07	7.992422e+04
3	bottom	4.254317e+08	-7.160426e+07	1.103976e+07
4	top	-1.016277e+08	-4.005265e+07	4.155850e+07
4	bottom	2.371977e+08	-7.007631e+07	3.682702e+07
5	top	-1.042836e+08	-4.963406e+07	-4.818647e+07
5	bottom	-1.045370e+08	-5.935929e+07	-5.473452e+07

Lamina	Location	eps1		eps2		gam12
1	top	-6.469278e-04	-3.637803e-03	-4.998701e-03
1	bottom	-6.330017e-04	-4.441136e-03	-6.486894e-03
2	top	-6.801947e-04	-4.393943e-03	6.541388e-03
2	bottom	-6.348065e-04	-5.228738e-03	7.993251e-03
3	top	1.677890e-03	-7.541434e-03	1.816459e-05
3	bottom	2.131908e-03	-8.784858e-03	-1.816459e-05
4	top	-5.894184e-04	-6.063532e-03	9.445114e-03
4	bottom	-5.440302e-04	-6.898327e-03	1.089698e-02
5	top	-5.912232e-04	-6.851134e-03	-1.095147e-02
5	bottom	-5.772971e-04	-7.654467e-03	-1.243966e-02

Smear Properties
Ex=	7.378204e+10	Ey=	1.346628e+10	Gxy=	2.604134e+10
vxy=	1.281868e+00	vyx=	2.339591e-01
vxz=	-7.111205e-01	vyz=	1.659425e-01
Nxyx=	1.296123e-16	Nxyy=	3.059588e-18	Nxxy=	-4.574660e-17	Nyxy=	5.916690e-18
ax=	-2.785747e-06	ay=	1.658942e-05	axy=	2.597093e-22	az=	1.301156e-05

Lamina	Location	fac1		fac2		fac12		Tsai Wu
1	top	8.958389e-02	2.795284e-01	-2.854233e-01	-4.686746e-01
1	bottom	8.958389e-02	2.795284e-01	-2.854233e-01	-4.686746e-01
2	top	9.170865e-02	3.274354e-01	3.517030e-01	-4.622865e-01
2	bottom	9.170865e-02	3.274354e-01	3.517030e-01	-4.622865e-01
3	top	2.037007e-01	5.019066e-01	-7.992422e-04	-2.650317e-01
3	bottom	2.037007e-01	5.019066e-01	-7.992422e-04	-2.650317e-01
4	top	8.442550e-02	4.275643e-01	4.794670e-01	-3.674531e-01
4	bottom	8.442550e-02	4.275643e-01	4.794670e-01	-3.674531e-01
5	top	9.039447e-02	4.740330e-01	-5.473452e-01	-2.800410e-01
5	bottom	9.039447e-02	4.740330e-01	-5.473452e-01	-2.800410e-01
```
