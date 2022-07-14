%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: 28/03/2022
% Author/s: Pol Padilla, Arturo Baltanas
%-------------------------------------------------------------------------%

clear;
close all;

%% INPUT DATA

% Material properties
E = 85e9; %[Pa]

% Cross-section parameters %[m]
t1 = 1.5e-3; 
t2 = 4e-3;
h1 = 500e-3;
h2 = 250e-3;
b = 775e-3;

% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35000;

% Number of elements for each part
nel = [3,6,12,24,48]; 
Nel = 96; % Number of elements for the "exact" solution

%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area
bt=sqrt(b^2+((h1-h2)/2)^2);   % Longitud barra inclinada
A=2*t1*bt+t2*h1+t2*h2;        % Area de la sección (no interior)
% Iz - Section inertia
centy=0;
theta=acos(b/bt);
centz=(t2*h1*b+2*t1*bt*b/2+t2*h2*0)/A;
centroide=[centz centy];

Izz1=1/12 *t1*bt^3*(sin(theta))^2;
Izz2=1/12 *t2*h1^3;
Izz3=1/12 *t1*bt^3*(sin(-theta))^2;
Izz4=1/12 *t2*h2^3;

Iyy1=1/12 *t1*bt^3*(cos(theta))^2;
Iyy2=0;
Iyy3=1/12 *t1*bt^3*(cos(theta))^2;
Iyy4=0;

Izy1=1/24 *t1*bt^3*sin(2*theta);
Izy2=0;
Izy3=1/24 *t1*bt^3*sin(-2*theta);
Izy4=0;

Iz=Izz1+Izz2+Izz3+Izz4+2*t1*bt*((h1+h2)/4)^2;
Iy=Iyy1+Iyy2+Iyy3+Iyy4+2*t1*bt*(b/2-centz)^2+t2*h1*(b-centz)^2;
Izy=Izy1+Izy2+Izy3+Izy4+t1*bt*((h1+h2)/4)*(b/2-centz)+t1*bt*(-(h1+h2)/4)*(b/2-centz);


% Compute parameter l:
% l - Equilibrium parameter
syms l x
eqLift1=l*(0.8-0.2*cos(pi*x/L1));
eqLift2=l*(1-(x-L1)/L2)*(1+(x-L1)/L2);
eqW1=M/(4*(L1+L2))+3*M/(2*L2^2) * (L1-x);
eqW2=M/(4*(L1+L2));
Lift=int(eqLift1,x,0,L1)+int(eqLift2,x,L1,L1+L2);
Weight=(int(eqW1,x,0,L1)+int(eqW2,x,L1,L1+L2)+Me)*g;
equil_l=vpasolve(Lift-Weight==0,l);

% Plot analytical solution
eqLift1=equil_l*(0.8-0.2*cos(pi*x/L1));
eqLift2=equil_l*(1-(x-L1)/L2)*(1+(x-L1)/L2);

    DeltaX=(L1+L2)/(100*Nel);
    x_an=0:DeltaX:L1+L2;

    Qy_vect=zeros(length(x_an),1);
    for i=1:find(x_an==L1)-1
        Qy_vect(i,1)=subs(eqLift1,x,x_an(i))-g*subs(eqW1,x,x_an(i));
    end
    for i=find(x_an==L1):length(x_an)
        Qy_vect(i,1)=subs(eqLift2,x,x_an(i))-g*subs(eqW2,x,x_an(i));
    end
    
    Fy_an=zeros(length(x_an),1);
    for i=length(x_an):-1:find(x_an==L1)
        Fy_an(i,1)=int(eqLift2,x,x_an(i),L1+L2)-g*int(eqW2,x,x_an(i),L1+L2);
    end
    for i=find(x_an==L1)-1:-1:1
        Fy_an(i,1)=int(eqLift1,x,x_an(i),L1)-g*int(eqW1,x,x_an(i),L1)+int(eqLift2,x,L1,L1+L2)-g*int(eqW2,x,L1,L1+L2)-g*Me;
    end

    %{
    Mz_an=zeros(length(x_an),1);
    for i=length(x_an):-1:find(x_an==L1)
        Mz_an(i,1)=int(eqLift2*x,x,x_an(i),L1+L2)-g*int(eqW2*x,x,x_an(i),L1+L2);
    end
    for i=find(x_an==L1)-1:-1:1
        Mz_an(i,1)=int(eqLift1*x,x,x_an(i),L1)-g*int(eqW1*x,x,x_an(i),L1)+int(eqLift2*x,x,L1,L1+L2)-g*int(eqW2*x,x,L1,L1+L2)-g*Me*(L1-x_an(i));
    end
    %}

    Mz_an=zeros(length(x_an),1);
    for i=length(x_an):-1:1
        Mz_an(i,1)=DeltaX*trapz(Fy_an(i:end,1));
    end
    
    u_an=zeros(length(x_an),1);
    theta_an=zeros(length(x_an),1);
   
    xB_x=zeros(length(x_an),length(x_an));
    for i=1:length(x_an)
    for j=1:length(x_an)
        xB_x(i,j)=x_an(i)-x_an(j);
    end
    end

    for i=1:length(x_an) % Teoremas Mohr
        theta_an(i,1)=DeltaX*trapz(Mz_an(1:i,1))/(E*Iz);
        u_an(i,1)=DeltaX*trapz(Mz_an(1:i,1).*transpose(xB_x(i,1:i)))/(E*Iz);
    end

fig = plotBeamsInitialize(L1+L2,x_an,u_an,theta_an,Fy_an,Mz_an);

uL1L2=zeros(length(nel),1);

%% Loop through each of the number of elements
for k = 1:length(nel)

    % PREPROCESS
    
    % Dimensions
    n_d = 1;                      % Number of dimensions
    n_el = nel(k);                % Total number of beams
    n_i = 2;                      % Number of DOFs for each node
    n_ne = 2;                     % Number of nodes in a beam
    n_nod = n_el+1;               % Total number of nodes
    n_dof = n_nod*n_i;            % Total number of degrees of freedom

    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];

    % Fix nodes matrix creation
    %  fixNod(k,1) = node at which some DOF is prescribed
    %  fixNod(k,2) = DOF prescribed (local, 1-2)
    %  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
    fixNod = [% Node        DOF  Magnitude
              % Write the data here...
              1 1 0
              1 2 0
    ];


    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    Tmat = zeros(n_el,1);
    Tmat(:,1)=1;

    % Nodal coordinates
    %  x(a,j) = coordinate of node a in the dimension j
    % Complete the coordinates
    x=zeros(n_nod,n_d);
    DeltaX=(L1+L2)/n_el;
    x(:,1)=transpose(0:DeltaX:L1+L2);

    % Nodal connectivities  
    %  Tnod(e,a) = global nodal number associated to node a of element e
    Tnod=zeros(n_el,n_ne);
    vectnod=[1 2];
    for i=1:n_el     
        Tnod(i,:)=vectnod;
        vectnod=vectnod+1;
    end

    % Computation of the DOFs connectivities
    Td = connectDOFs(n_el,n_ne,n_i,Tnod);

    % Computation of element stiffness matrices
    Kel = computeKelBar(n_ne,n_i,n_el,x,Tnod,mat,Tmat);

    % Computation element force vector
    F_el = computeF(n_el,n_ne,n_i,Qy_vect,x_an,Tmat,mat,x,Tnod,Me,g,L1);

    % Global matrix assembly & force vector assembly
    [KG,Fext] = assembly_KG_Fext(n_el,n_dof,n_ne,n_i,Td,Kel,F_el);

    % Apply conditions 
    [vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

    % System resolution
    [u,R] = solveSys(vL,vR,uR,KG,Fext);

    % Compute internal distributions
    [pu,pt,Fy,Mz] = computeInternal(n_el,n_ne,n_i,x,u,Kel,Td,Tnod,mat,Tmat);

    
    %% SOLVER
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]
    
    %% POSTPROCESS
    
    uL1L2(k,1)=u(n_dof-1,1);

    % Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
end

% Add figure legends
sLeg=zeros(1,length(nel)+1);
    sLeg(1,1)=Nel;
    sLeg(1,2:end)=nel;

figure(fig)
legend(strcat('N=',cellstr(string(sLeg))),'location','northeast');

% Convergence error
Er=abs((u_an(Nel+1,1)-uL1L2)./uL1L2);
    figure (2);
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    semilogx(nel,Er)
    title('Convergence error, displacements at wing tip ')
    xlabel("Number of elements")
    ylabel("Error $\varepsilon_r$")
    grid on
    grid minor

%% Von Mises Criterion
% -------------------------------------------------------------------------
% Analisis for all sections and all points in the section. Consumes a lot 
% of time, stop executing after a minute to see the results of the first 
% sections.
% -------------------------------------------------------------------------

VonMisses=zeros(length(x_an),1);
position_s_local_max=zeros(length(x_an),3); % [se12 se23 se34]
localization_max=zeros(length(x_an),1); % if it's section 12, 23, 34

for i=1:length(x_an)
y=h1/2; x=b-centz;
syms y
sigma_x=Iy*y*Mz_an(i,1)/(Iz*Iy);

syms s
q12o= -Fy_an(i,1)/Iz * s^2/2 *t2;
    q2=-Fy_an(i,1)/Iz * h1^2/8 *t2;
q23o= -Fy_an(i,1)/Iz *s*t1* (h1/2 -0.5*s*sin(theta)) +q2;
    q3=-Fy_an(i,1)/Iz *bt*t1* (h1/2 -0.5*bt*sin(theta)) +q2;
q34o= -Fy_an(i,1)/Iz *s*t2* (h2/2-s/2) +q3;

r12=b-centz;
r23=(tan(theta)*centz+h2/2)/sqrt(1+(tan(theta))^2);
r34=centz;

intq12=int(q12o*r12,s,0,h1/2);
intq23=int(q23o*r23,s,0,bt);
intq34=int(q34o*r34,s,0,h2/2);
sumint=intq12+intq23+intq34;

Ain=(h1+h2)/2 * b;
qs0=2*sumint/(2*Ain);

q12= q12o+qs0;    
q23= q23o+qs0;   
q34= q34o+qs0;

% seccion 12
se=0:0.001:h1/2;
VM12=zeros(length(se),1);
for c=1:length(se)
    tau12=subs(q12,s,se(c))/t2;
    ye=se(c);
    sigma12=subs(sigma_x,y,ye);
    VM12(c,1)=sqrt(sigma12^2+3*tau12^2);
end
position_s_local_max(i,1)=se(find(VM12==max(VM12)));

% seccion 23
se=0:0.001:bt;
VM23=zeros(length(se),1);
for c=1:length(se)
    tau23=subs(q23,s,se(c))/t1;
    ye=h1/2-sin(theta)*se(c);
    sigma23=subs(sigma_x,y,ye);
    VM23(c,1)=sqrt(sigma23^2+3*tau23^2);
end
position_s_local_max(i,2)=se(find(VM23==max(VM23)));

% seccion 34
se=0:0.001:h2/2;
VM34=zeros(length(se),1);
for c=1:length(se)
    tau34=subs(q34,s,se(c))/t2;
    ye=h2/2-se(c);
    sigma34=subs(sigma_x,y,ye);
    VM34(c,1)=sqrt(sigma34^2+3*tau34^2);
end
position_s_local_max(i,3)=se(find(VM34==max(VM34)));

VonMisses(i,1)=max([max(VM12) max(VM23) max(VM34)]);
    if VonMisses(i,1)==max(VM12)
        localization_max(i,1)=12;
    end
    if VonMisses(i,1)==max(VM23)
        localization_max(i,1)=23;
    end
    if VonMisses(i,1)==max(VM34)
        localization_max(i,1)=34;
    end

end



figure
hold on
title('Von misses a lo largo del ala, 'interpreter', 'latex', 'FontSize',16)
plot(x_an,VonMisses, 'Color', orange, 'LineWidth',1.2)
xlabel('x[m]', 'Interpreter','latex', 'FontSize',15)
ylabel('$\sigma$[MPa]', 'Interpreter','latex', 'FontSize',15)
grid on
grid minor
axis padded
hold off


