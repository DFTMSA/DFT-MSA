%==========================================================================
%DFT-MSA solver
%==========================================================================
%AUTHOR
%P. Cornelissen

%MODEL DESCRIPTION
%This function solves the Poisson problem including ion-ion correlation
%between two charged plates, using density %functional theory (DFT) and the
%mean sphere approximation (MSA). This allows computation of the ion
%concentration profiles and electric potential in the nanopore. This model
%was used for the study of Cornelissen et al. (2021). The numerical 
%procedure is based on Le et al. (2015) and Rocha et al. (2016). 

%INPUT PARAMETERS
%ncell = number of nodes for the one-dimensional finite difference discretization of the
%pore space
%Cb_Cl = chloride concentration [Cl-] in the bulk solution [mol/m^3]
%Cb_Na = sodium concentration [Na+] in the bulk solution [mol/m^3]
%sigma = surface charge density of the clay mineral [C/m^2]
%d = ion diameter [m]
%h = pore size [m]

%OUTPUT PARAMETERS
%x = x-coordinate [m]
%C_Cl = chloride concentration in pore [mol/m^3]
%C_Na = sodium concentration in pore [mol/m^3]
%C_Ca = calcium concentration in pore [mol/m^3]
%g_Cl = chloride ion-ion correlation term in pore [-]
%g_Na = sodium ion-ion correlation term in pore [-]
%g_Ca = calcium ion-ion correlation term in pore [-]
%psi = electric potential in pore [V]

%REFERENCES
%Cornelissen, P., van der Zee, S.E.A.T.M., Leijnse, A., Niasar, V., (2021).
   %Submitted manuscript

%Le, T.D., Rocha, A.C., Murad, M.A., Moyne, C., Oliveira, S.P., (2016). Multiscale contaminant transport in swelling clays incorporating ion-ion 
   %correlation effects. In: Handbook of Groundwater Engineering, 3rd ed. Editors: Cushman, J.H., Tartakovsky, D.M. CRC Press, Taylor & Francis
   %Group, Boca Raton. https://doi.org/10.1201/9781315371801

%Rocha, A.C., Murad, M.A., Moyne, C., Oliveira, S.P., Le, T.D., (2016). A new methodology for computing ionic profiles and disjoining pressure 
   %in swelling porous media. Computional Geosciences, 20, 975-996. https://doi.org/10.1007/s10596-016-9572-5 

function [x,C_Cl,C_Na,C_Ca,g_Cl,g_Na,g_Ca,psi] = DFTMSA_solver(ncell,Cb_Cl,Cb_Na,sigma,d,h)
if Cb_Na > Cb_Cl
    msg = 'Error: Bulk sodium concentration should be smaller than the bulk chloride concentration!';
    error(msg);
end
if h < 0.5*d
    msg = 'Error: pore size should be larger than the ion radius!';
    error(msg);
end

%Physical constants
F = 96485.33289; %Faraday constant [C/mol]
R = 8.3144598; %Gas constant [J/mol/K]
nA = 6.02214076e23; %Avogadro's number [1/mol]
e = F/nA; %Electron charge [C]
kB = R/nA; %Boltzmann constant [J/K]
perm0 = 8.854e-12; %Vacuum permittivity [C^2/J/m]

%Variables
Cb_Ca = 0.5*(Cb_Cl - Cb_Na); %[Ca2+] = calcium concentration in bulk solution [mol/m^3]
T = 298; %Temperature [K]
perm_r = 78.5; %Relative permittivity water [-]
perm = perm0*perm_r; %Absolute permittivity water [C^2/J/m]
beta = 1/(kB*T);

%Solver options
tol = 1e-6; %Tolerance for linear solver
maxit = 20; %Maximum iterations for linear solver

%Grid initialization
dx = (h-d/2)/(ncell-1); %Distance between nodes [m]
%Calculate x-coordinate of each node
x = zeros(ncell,1);
for i = 2:ncell
    x(i) = x(i-1) + dx;
end

%First, the coefficients containing the ion-ion correlation are calculated,
%following Appendix A.3. of Cornelissen et al. (2021).
%Initialization of the kernels
K1 = zeros(ncell,ncell); %Kernel describing Na-Na and Cl-Cl interaction
K2 = zeros(ncell,ncell); %Kernel describing Na-Cl and Cl-Na interaction
K3 = zeros(ncell,ncell); %Kernel describing Cl-Ca and Ca-Cl interaction
K4 = zeros(ncell,ncell); %Kernel describing Na-Ca and Ca-Na interaction
K5 = zeros(ncell,ncell); %Kernel describing Ca-Ca interaction

%Initialization of the ion-ion correlation terms
M1 = zeros(ncell,ncell); %Hard-sphere component for x' > 0
M2 = zeros(ncell,ncell); %Hard-sphere component for x'< 0
L1 = zeros(ncell,ncell); %Short-range Coulomb correction for x' > 0
L2 = zeros(ncell,ncell); %Short-range Coulomb correction for x' < 0
J = zeros(ncell,1); %Hard sphere component within the exclusion zone

rho = nA*(Cb_Cl + Cb_Na + Cb_Ca); %Total ion concentration in bulk solution
eta = rho*pi*d^3/6; %Packing factor
a1 = -(1+2*eta)^2/(1-eta)^4; %Coefficient for hard-sphere component
a2 = 6*eta*(1+0.5*eta)^2/(1-eta)^4; %Coefficient for hard-sphere component
a3 = eta*a1/2; %Coefficient for hard-sphere component
y = sqrt(F^2*(Cb_Cl+Cb_Na+4*Cb_Ca)/(perm*R*T))*d; %Inverse Debye length
B = (y^2+y-y*sqrt(1+2*y))/(y^2); %Coefficient for the short-range Coulomb correction

a = h - 0.5*d; %x-coordinate of start of the exclusion zone
b = h + 0.5*d; %x-coordinate of end of the exclusion zone

%Computation of the kernels
for i = 1:ncell
    for j = 1:ncell
        r1 = x(i)+x(j); %Distance |x-x'| for x' < 0
        r2 = abs(x(i)-x(j)); %Distance |x-x'| for x' > 0        
        %If x' < 0
        if x(i) < d && (x(i)+x(j)) < d
            M1(i,j) = 2*pi*(a1*(d^2/2-r1^2/2)+a2/d*(d^3/3-r1^3/3) + a3/d^3*(d^5/5-r1^5/5));%Hard-sphere component for x' > 0
            L1(i,j) = d*(d-r1) - B*(d^2-r1^2) + B^2/(3*d)*(d^3-r1^3); %Short-range Coulomb correction for x' > 0
            K1(i,j) = K1(i,j) + nA*(M1(i,j) + beta*e^2/(2*perm*d)*L1(i,j)); %Kernel describing Na-Na and Cl-Cl interaction
            K2(i,j) = K2(i,j) + nA*(M1(i,j) - beta*e^2/(2*perm*d)*L1(i,j)); %Kernel describing Na-Cl and Cl-Na interaction
            K3(i,j) = K3(i,j) + nA*(M1(i,j) - beta*e^2/(perm*d)*L1(i,j)); %Kernel describing Cl-Ca and Ca-Cl interaction
            K4(i,j) = K4(i,j) + nA*(M1(i,j) + beta*e^2/(perm*d)*L1(i,j)); %Kernel describing Na-Ca and Ca-Na interaction
            K5(i,j) = K5(i,j) + nA*(M1(i,j) + 2*beta*e^2/(perm*d)*L1(i,j)); %Kernel describing Ca-Ca interaction
        end
        %If x' > 0 
        if (x(i)-d) < x(j) && x(j) < (x(i)+d)
            M2(i,j) = 2*pi*(a1*(d^2/2-r2^2/2)+a2/d*(d^3/3-r2^3/3) + a3/d^3*(d^5/5-r2^5/5));%Hard-sphere component for x'< 0
            L2(i,j) = d*(d-r2) - B*(d^2-r2^2) + B^2/(3*d)*(d^3-r2^3); %Short-range Coulomb correction for x' < 0            
            K1(i,j) = K1(i,j) + nA*(M2(i,j) + beta*e^2/(2*perm*d)*L2(i,j)); %Kernel describing Na-Na and Cl-Cl interaction
            K2(i,j) = K2(i,j) + nA*(M2(i,j) - beta*e^2/(2*perm*d)*L2(i,j)); %Kernel describing Na-Cl and Cl-Na interaction
            K3(i,j) = K3(i,j) + nA*(M2(i,j) - beta*e^2/(perm*d)*L2(i,j)); %Kernel describing Cl-Ca and Ca-Cl interaction
            K4(i,j) = K4(i,j) + nA*(M2(i,j) + beta*e^2/(perm*d)*L2(i,j)); %Kernel describing Na-Ca and Ca-Na interaction
            K5(i,j) = K5(i,j) + nA*(M2(i,j) + 2*beta*e^2/(perm*d)*L2(i,j)); %Kernel describing Ca-Ca interaction
        end
    end
    %Computation of the hard sphere component within the exclusion zone
    if h > 0.5*d && h < 1.5*d && x(i) < min(1.5*d-h,h-0.5*d)
        tmp1 = 0.5*a1*(d^2*(b-a)-(1/3)*((x(i)+b)^3-(x(i)+a)^3));
        tmp2 = a2/(3*d)*(d^3*(b-a)-0.25*((x(i)+b)^4-(x(i)+a)^4));
        tmp3 = a3/(5*d^3)*(d^5*(b-a)-(1/6)*((x(i)+b)^6-(x(i)+a)^6));
        J(i) = J(i) + tmp1 + tmp2 + tmp3; %Combine
    end
    if h > 0.5*d && x(i) > max(0,h-1.5*d)  && x(i) < h - 0.5*d
        tmp1 = 0.5*a1*(d^2*(b-a)-(x(i)^2*(b-a)-x(i)*(b^2-a^2)+(1/3)*(b^3-a^3)));
        tmp2 = a2/(3*d)*(d^3*(b-a)+0.25*((x(i)-b)^4-(x(i)-a)^4));
        tmp3 = a3/(5*d^3)*(d^5*(b-a)+(1/6)*((x(i)-b)^6-(x(i)-a)^6));
        J(i) = J(i) + tmp1 + tmp2 + tmp3;
    end
    J(i) = nA*(Cb_Cl + Cb_Na + Cb_Ca)*J(i); %Multiply with bulk concentrations
end

%Compute the weights obtained by discretizing the integral equation with
%the trapezoidal rule
w = dx*ones(ncell,1);
w(1) = dx/2;
w(ncell) = dx/2;

%Setup the eigenvalue problem
A1 = zeros(ncell,ncell);
A2 = zeros(ncell,ncell);
A3 = zeros(ncell,ncell);
A4 = zeros(ncell,ncell);
A5 = zeros(ncell,ncell);
for i = 1:ncell
    for j = 1:ncell
        A1(i,j) = w(j)*K1(i,j);
        A2(i,j) = w(j)*K2(i,j);
        A3(i,j) = w(j)*K3(i,j);
        A4(i,j) = w(j)*K4(i,j);
        A5(i,j) = w(j)*K5(i,j);
    end
end

%Solve the eigenvalue problem (Equation A9 in Cornelissen et al., 2021)
[U1,lambda1] = eig(A1); %Eigenvalue problem for Na-Na and Cl-Cl interaction
[U2,lambda2] = eig(A2); %Eigenvalue problem for Na-Cl and Na-Cl interaction
[U3,lambda3] = eig(A3); %Eigenvalue problem for Na-Ca and Ca-Na interaction
[U4,lambda4] = eig(A4); %Eigenvalue problem for Cl-Ca and Ca-Cl interaction
[U5,lambda5] = eig(A5); %Eigenvalue problem for Ca-Ca interaction
%Convert the eigenvalue-matrices to vectors
lambda1 = diag(lambda1); %Eigenvectors for Na-Na and Cl-Cl interaction
lambda2 = diag(lambda2);
lambda3 = diag(lambda3);
lambda4 = diag(lambda4);
lambda5 = diag(lambda5);

%Compute the rotation matrices
I_12 = inv(U2)*U1;
I_13 = inv(U3)*U1;
I_14 = inv(U4)*U1;
I_15 = inv(U5)*U1;

% First iteration: Poisson-Boltzmann equation
%For the first iteration, we neglect ion-ion correlations, and thus solve
%the nonlinear Poisson-Boltzmann equation
%The nonlinear Poisson-Boltzmann equation is linearize using Newton's
%method. This results in a linear system of equations Ax = b, which is
%solved for x (the electric potential) in an iterative fashion.

%1a. Calculate 
%Analytical solution for the linearized Poisson-Boltzmann equation
k = sqrt((Cb_Cl + Cb_Na + 4*Cb_Ca)*F^2/(perm*R*T)); %Inverse Debye length for analytical solution
psi_analytical = sigma/(perm*k)*cosh(k*x)/sinh(k*(h-0.5*d)); %Analytical solution for the electric potential

%Initialize parameters
psi_old = psi_analytical; %Initial guess
it = 0; %Iteration counter
delta_psi = 1; %Change in solution
k = sqrt(F^2/(perm*R*T)); %Inverse Debye length for numerical solution
A_matrix = zeros(ncell,ncell); %A matrix in linear system Ax = b
b_matrix = zeros(ncell,1); %b matrix in linear system Ax = b
%Start solver
while delta_psi > 1e-6
    %Step 1: Setup of A and b matrices
    %Boundary conditions (outer nodes)
    %x = 0: d/dx(psi) = 0
    A_matrix(1,1) = -(2 + Cb_Cl*k^2*dx^2*exp(psi_old(1)) + Cb_Na*k^2*dx^2*exp(-psi_old(1)) + 4*Cb_Ca*k^2*dx^2*exp(-2*psi_old(1)));
    A_matrix(1,2) = 2;
    b_matrix(1) = Cb_Cl*k^2*dx^2*exp(psi_old(1))*(1 - psi_old(1)) - Cb_Na*k^2*dx^2*exp(-psi_old(1))*(1 + psi_old(1)) - 2*Cb_Ca*k^2*dx^2*exp(-2*psi_old(1))*(1 + 2*psi_old(1));
    %x = h - 0.5*d: d/dx(psi) = -sigma*F/(perm*R*T)
    A_matrix(ncell,ncell) = -(2 + Cb_Cl*k^2*dx^2*exp(psi_old(ncell)) + Cb_Na*k^2*dx^2*exp(-psi_old(ncell)) + 4*Cb_Ca*k^2*dx^2*exp(-2*psi_old(ncell)));
    A_matrix(ncell,ncell-1) = 2;
    b_matrix(ncell) = Cb_Cl*k^2*dx^2*exp(psi_old(ncell))*(1 - psi_old(ncell)) - Cb_Na*k^2*dx^2*exp(-psi_old(ncell))*(1 + psi_old(ncell)) - 2*Cb_Ca*k^2*dx^2*exp(-2*psi_old(ncell))*(1 + 2*psi_old(ncell)) - 2*F*dx*sigma/(perm*R*T);
    %Inner nodes: Poisson equation discretized with finite difference
    %scheme
    for i = 2:(ncell-1)
        A_matrix(i,i) = -(2 + Cb_Cl*k^2*dx^2*exp(psi_old(i)) + Cb_Na*k^2*dx^2*exp(-psi_old(i)) + 4*Cb_Ca*k^2*dx^2*exp(-2*psi_old(i)));
        A_matrix(i,i-1) = 1;
        A_matrix(i,i+1) = 1;
        b_matrix(i) = Cb_Cl*k^2*dx^2*exp(psi_old(i))*(1 - psi_old(i)) - Cb_Na*k^2*dx^2*exp(-psi_old(i))*(1 + psi_old(i)) - 2*Cb_Ca*k^2*dx^2*exp(-2*psi_old(i))*(1 + 2*psi_old(i));
    end
    %Step 2: Solve system of equations
    %LU-decomposition of the A-matrix
    [L,U] = lu(A_matrix);
    %Solve linear system of equations to obtain electric potential
    [psi,flag] = bicgstab(A_matrix,b_matrix,tol,maxit,L,U);
    %Step 3: Check convergence criteria
    %Compute change in solution
    delta_psi = sqrt(sum((psi-psi_old).^2));
    %Update electric potential for new iteration
    psi_old = psi;
end
%The result is the electric potential as calculated by Poisson-Boltzmann
%theory
psi_pb = psi;

%Second iteration onwards: include ion-ion correlation
%Initializaton of variables
a1 = ones(ncell,1);
b1 = ones(ncell,1);
c1 = ones(ncell,1);
delta_psi2 = 1;
psi_k = psi;

while delta_psi2 > 1e-6 %Iterate over k for convergence in staggered algorithm
    psi_k = psi;
    
    %STEP 1: calculate the coefficients a1, b1, and c1 (Equation A11 in
    %Cornelissen et al., 2021). The nonlinear system of equations is
    %linearized using a Picard-style iteration algorithm, resuling in a
    %linear system of equations to be solved (Ax = b).
    %Update coefficients from previous iteration
    a1_old = a1;
    b1_old = b1;
    c1_old = c1;
    
    %Initialize A and b matrices
    A_matrix = zeros(3*ncell,3*ncell);
    b_matrix = zeros(3*ncell,1);
    
    %Compute A and b matrices
    A_matrix(1:ncell,1:ncell) = U1;
    A_matrix(ncell+1:2*ncell,ncell+1:2*ncell) = U1;
    A_matrix(2*ncell+1:3*ncell,2*ncell+1:3*ncell) = U1;
    b_matrix(1:ncell) = exp(psi + Cb_Cl*U1*(a1_old.*lambda1) + Cb_Na*U2*((I_12*b1_old).*lambda2) + Cb_Ca*U3*((I_13*c1_old).*lambda3) + J) - 1;
    b_matrix(ncell+1:2*ncell) = exp(-psi + Cb_Cl*U2*((I_12*a1_old).*lambda2) + Cb_Na*U1*(b1_old.*lambda1) + Cb_Ca*U4*((I_14*c1_old).*lambda4) + J) - 1;
    b_matrix(2*ncell+1:3*ncell) = exp(-2*psi + Cb_Cl*U3*((I_13*a1_old).*lambda3) + Cb_Na*U4*((I_14*b1_old).*lambda4) + Cb_Ca*U5*((I_15*c1_old).*lambda5) + J) - 1; 
    %Solve linear system of equations
    %LU-decomposition of the A-matrix
    [L,U] = lu(A_matrix);
    %Solve linear system of equations to obtain coefficients
    [out,flag] = bicgstab(A_matrix,b_matrix,tol,maxit,L,U);    
    %Extract a1, b1, and c1 for new iteration from solution
    a1 = out(1:ncell);
    b1 = out(ncell+1:2*ncell);
    c1 = out(2*ncell+1:3*ncell);
    
    %STEP 2: compute the ion-ion correlation term in the ion distribution
    %from new coefficients a1, b1, and c1
    g_Cl = exp(Cb_Cl*U1*(a1.*lambda1) + Cb_Na*U2*((I_12*b1).*lambda2) + Cb_Ca*U3*((I_13*c1).*lambda3) + J);
    g_Na = exp(Cb_Cl*U2*((I_12*a1).*lambda2) + Cb_Na*U1*(b1.*lambda1) + Cb_Ca*U4*((I_14*c1).*lambda4) + J);
    g_Ca = exp(Cb_Cl*U3*((I_13*a1).*lambda3) + Cb_Na*U4*((I_14*b1).*lambda4) + Cb_Ca*U5*((I_15*c1).*lambda5) + J);
    
    %STEP 3: compute the electric potential using the updated ion-ion
    %correlation terms
    psi_old = psi;
    delta_psi = 1;
    
    while delta_psi > 1e-6 %Iterate over s convergence in Newton scheme for Poisson problem   
        %Step 1: Setup of A and b matrices
        %Boundary conditions (outer nodes)
        %x = 0: d/dx(psi) = 0
        A_matrix = zeros(ncell,ncell);
        b_matrix = zeros(ncell,1);
        
        A_matrix(1,1) = -(2 + Cb_Cl*k^2*dx^2*g_Cl(1)*exp(psi_old(1)) + Cb_Na*k^2*dx^2*g_Na(1)*exp(-psi_old(1)) + 4*Cb_Ca*k^2*dx^2*g_Ca(1)*exp(-2*psi_old(1)));
        A_matrix(1,2) = 2;
        b_matrix(1) = Cb_Cl*k^2*dx^2*g_Cl(1)*exp(psi_old(1))*(1 - psi_old(1)) - Cb_Na*k^2*dx^2*g_Na(1)*exp(-psi_old(1))*(1 + psi_old(1)) - 2*Cb_Ca*k^2*dx^2*g_Ca(1)*exp(-2*psi_old(1))*(1 + 2*psi_old(1));
        
        A_matrix(ncell,ncell) = -(2 + Cb_Cl*k^2*dx^2*g_Cl(ncell)*exp(psi_old(ncell)) + Cb_Na*k^2*dx^2*g_Na(ncell)*exp(-psi_old(ncell)) + 4*Cb_Ca*k^2*dx^2*g_Ca(ncell)*exp(-2*psi_old(ncell)));
        A_matrix(ncell,ncell-1) = 2;
        b_matrix(ncell) = Cb_Cl*k^2*dx^2*g_Cl(ncell)*exp(psi_old(ncell))*(1 - psi_old(ncell)) - Cb_Na*k^2*dx^2*g_Na(ncell)*exp(-psi_old(ncell))*(1 + psi_old(ncell)) - 2*Cb_Ca*k^2*dx^2*g_Ca(ncell)*exp(-2*psi_old(ncell))*(1 + 2*psi_old(ncell)) - 2*F*dx*sigma/(perm*R*T);
        for i = 2:(ncell-1)
            A_matrix(i,i) = -(2 + Cb_Cl*k^2*dx^2*g_Cl(i)*exp(psi_old(i)) + Cb_Na*k^2*dx^2*g_Na(i)*exp(-psi_old(i)) + 4*Cb_Ca*k^2*dx^2*g_Ca(i)*exp(-2*psi_old(i)));
            A_matrix(i,i-1) = 1;
            A_matrix(i,i+1) = 1;
            b_matrix(i) = Cb_Cl*k^2*dx^2*g_Cl(i)*exp(psi_old(i))*(1 - psi_old(i)) - Cb_Na*k^2*dx^2*g_Na(i)*exp(-psi_old(i))*(1 + psi_old(i)) - 2*Cb_Ca*k^2*dx^2*g_Ca(i)*exp(-2*psi_old(i))*(1 + 2*psi_old(i));
        end
        
        %Solve linear system of equations
        %LU-decomposition of the A-matrix
        [L,U] = lu(A_matrix);
        %Solve linear system of equations to obtain electric potential
        [psi,flag] = bicgstab(A_matrix,b_matrix,tol,maxit,L,U);
        %Compute change in solution compared to previous iteration in the
        %Newton scheme
        delta_psi = sqrt(sum((psi-psi_old).^2));
        psi_old = psi;
    end   
    %Compute change in solution compared to previous iteration in the
    %staggered scheme
    delta_psi2 = sqrt(sum((psi-psi_k).^2));
end

%Now electric potential and ion-ion correlation terms are known
%Compute the ion concentration profiles
C_Cl = Cb_Cl*exp(psi).*g_Cl;
C_Na = Cb_Na*exp(-psi).*g_Na;
C_Ca = Cb_Ca*exp(-2*psi).*g_Ca;

%Convert electric potential units to V (volts)
psi = R*T*psi/F;

end
