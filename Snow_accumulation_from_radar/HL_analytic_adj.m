function [z, rho] = HL_analytic_adj(T, acc, rho_0)
%T = (273+(-40));						% firn temperature (Kelvin)
%acc = 0.12;					% ICE-equiv accumulation rate (m/year)
SEP = 0;						% distance Tx to Rx (metres)
%----------------------------------------------------------------------------
% (2) constants
% velocities c (metre per second) and relative permittivity E
c = 3e8; c_ice = 1.685e8; E_ice = (c/c_ice)^2;
% ice density (Mg/m^3)
rho_i = 0.917;
% firn density at surface
%rho_0 = 0.35;
% critical density (stage 1-2 densification transition)
rho_c = 0.55;
% water density
rho_w = 1;
%----------------------------------------------------------------------------
% (3) density  vs. depth profile (look-up table to 5000 metres)
% from Herron-Langway model
% rho vs. z, where rho is in Mg per cubic metres
z = [0:0.05:5000]';
R = 8.314;
k0 = 11*exp(-10160/(R*T));
k1 = 575*exp(-21400/(R*T));
alpha = log(rho_0/(rho_i-rho_0));
beta = log(rho_c/(rho_i-rho_c));
% WATER-equiv accumulation rate (m/year)
A = acc*rho_i/rho_w;
% critical depth hc (metre)
hc = (beta-alpha)/(rho_i*k0);
% density profile evaluation
Z=(z<hc).*exp(rho_i*k0*z+alpha) + (z>=hc).*exp(rho_i*k1*(z-hc)/(A^0.5)+beta);
rho=rho_i*Z./(1+Z);
end