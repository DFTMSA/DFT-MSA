%Example script for use of DFTMSA_solver.m
%Input parameters
ncell = 1001; %Define the number of cells for the finite difference discretization
Cb_Cl = 10; %[Cl-] = chloride concentration in bulk solution [mol/m^3]
Cb_Na = 5; %[Na+] = sodium concentration in bulk solution [mol/m^3]
sigma = -0.1; %Surface charge density of the clay mineral [C/m^2]
d = 4.25e-10; %Ion diameter [-]
h = d; %Interlameller half-spacing [m]

%Run the DFT-MSA model
[x,C_Cl,C_Na,C_Ca,g_Cl,g_Na,g_Ca,psi] = DFTMSA_solver(ncell,Cb_Cl,Cb_Na,sigma,d,h);

%Plot output
clf;
subplot(2,2,1)
plot(x/d,C_Cl)
xlabel('Distance [-]')
ylabel('Chloride concentration [mmol/l]')

subplot(2,2,2)
plot(x/d,C_Na)
xlabel('Distance [-]')
ylabel('Sodium concentration [mmol/l]')

subplot(2,2,3)
plot(x/d,C_Ca)
xlabel('Distance [-]')
ylabel('Calcium concentration [mmol/l]')

subplot(2,2,4)
plot(x/d,psi)
xlabel('Distance [-]')
ylabel('Electric potential [V]')
