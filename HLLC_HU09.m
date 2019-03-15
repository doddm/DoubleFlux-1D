function [ f_HLLC ] = HLLC_HU09( V_L,V_R,gammaS,e0S,side )
% HLLC_HU09 Compute flux with HLLC Riemann solver for interface interaction
% in compressible multi-fluid flow [X. Y. Hu, N. A. Adams, G. Iaccarino JCP 2009]
sign = @(s) (s>=0) - (s<=0);
      
u_L = V_L(2,:);
u_R = V_R(2,:);

rho_L = V_L(1,:);
rho_R = V_R(1,:);

p_L = V_L(3,:);
p_R = V_R(3,:);

a_L = sqrt(gammaS.*p_L./rho_L);
a_R = sqrt(gammaS.*p_R./rho_R);

T_L = getTfromPandRho( p_L,rho_L );  % temperature [K]
T_R = getTfromPandRho( p_R,rho_R );

E_L = p_L./(gammaS-1) + rho_L.*e0S + 0.5*rho_L.*u_L.^2; % Eq. (21)
E_R = p_R./(gammaS-1) + rho_R.*e0S + 0.5*rho_R.*u_R.^2; % Eq. (21)

e_L = getEnergyfromTandRho( T_L,rho_L ); % specific internal energy [J/kg]
e_R = getEnergyfromTandRho( T_R,rho_R );

u_roe = (sqrt(rho_L) .* u_L + sqrt(rho_R) .* u_R)./(sqrt(rho_L) + sqrt(rho_R));
% e_roe = (sqrt(rho_L) .* e_L + sqrt(rho_R) .* e_R)./(sqrt(rho_L) + sqrt(rho_R));
rho_roe = sqrt(rho_L .* rho_R);

p_rho_roe = (sqrt(rho_L) .* (p_L./rho_L) + sqrt(rho_R) .* (p_R./rho_R))./(sqrt(rho_L) + sqrt(rho_R)) + 0.5*((u_R - u_L) ./(sqrt(rho_L) + sqrt(rho_R))).^2;  % Eq. (18) X. Y. Hu et al. JCP 2009

H_L = p_L./rho_L + e_L + 0.5*u_L.^2;  % enthalpy [J/kg]
H_R = p_R./rho_R + e_R + 0.5*u_R.^2;

% H_roe = (sqrt(rho_L) .* H_L + sqrt(rho_R) .* H_R)./(sqrt(rho_L) + sqrt(rho_R));

% p_i = (gammaS-1).*rho_L;
% p_rho   = (gammaS-1).*(e_roe-e0S);

psi_L = (e_L-e0S).*(gammaS - 1);
psi_R = (e_R-e0S).*(gammaS - 1);
gamma = gammaS - 1;
gamma_roe = gamma;
psi_roe = (sqrt(rho_L) .* psi_L + sqrt(rho_R) .* psi_R)./(sqrt(rho_L) + sqrt(rho_R));

c_roe = sqrt(psi_roe + gamma_roe.*p_rho_roe);  % Eq. (16) X. Y. Hu et al. JCP 2009

% Calculate ghost states
delp = p_R - p_L;
delrho = rho_R - rho_L;

erg = e_L + (delp - psi_roe.*delrho)./(gamma_roe .* rho_roe); % Eq. (26) X. Y. Hu et al. JCP 2009
elg = e_R - (delp - psi_roe.*delrho)./(gamma_roe .* rho_roe);

% E_L = p_L./(gammaS-1) + rho_L.*e0S + 0.5*rho_L.*u_L.^2; % Eq. (21)
% E_R = p_R./(gammaS-1) + rho_R.*e0S + 0.5*rho_R.*u_R.^2; % Eq. (21)

if(side==0)
    E_L = rho_L.*elg + 0.5*rho_L.*u_L.^2;
    E_R = rho_R.*e_R + 0.5*rho_R.*u_R.^2; 
end
if(side==1)
    E_L = rho_L.*e_L + 0.5*rho_L.*u_L.^2; 
    E_R = rho_R.*erg + 0.5*rho_R.*u_R.^2;
end

% E_L = rho_L.*e_L + 0.5*rho_L.*u_L.^2; 
% E_R = rho_R.*e_R + 0.5*rho_R.*u_R.^2;

% E_L = rho_L.*elg + 0.5*rho_L.*u_L.^2;
% E_R = rho_R.*erg + 0.5*rho_R.*u_R.^2; 

% E_L = p_L./(gammaS-1) + rho_L.*e0S + 0.5*rho_L.*u_L.^2; % Eq. (21)
% E_R = p_R./(gammaS-1) + rho_R.*e0S + 0.5*rho_R.*u_R.^2; % Eq. (21)

f_L = [rho_L.*u_L; rho_L.*(u_L).^2+p_L; u_L.*(E_L+p_L)];
f_R = [rho_R.*u_R; rho_R.*(u_R).^2+p_R; u_R.*(E_R+p_R)];

W_L = [rho_L; rho_L.*u_L; E_L];
W_R = [rho_R; rho_R.*u_R; E_R];

s_L = min(u_L - a_L, u_roe - c_roe);  % Eq. (12) X. Y. Hu et al. JCP 2009
s_R = max(u_R + a_R, u_roe + c_roe);

% s_L = u_L - a_L;
% s_R = u_R + a_R;

sM = (p_L - p_R - rho_L.*u_L.*(s_L-u_L) + rho_R.*u_R.*(s_R-u_R))./(rho_R.*(s_R-u_R)-rho_L.*(s_L-u_L));
pM = rho_R.*(u_R-s_R).*(u_R-sM) + p_R;

WML = [(s_L-u_L)./(s_L-sM).*rho_L;
       ((s_L-u_L).*rho_L.*u_L + (pM-p_L))./(s_L-sM);
       ((s_L-u_L).*E_L - p_L.*u_L + pM.*sM)./(s_L-sM)];

WMR = [(s_R-u_R)./(s_R-sM).*rho_R;
       ((s_R-u_R).*rho_R.*u_R + (pM-p_R))./(s_R-sM);
       ((s_R-u_R).*E_R - p_R.*u_R + pM.*sM)./(s_R-sM)];

s_minus = min(0,s_L);
s_plus = max(0,s_R);

f_HLLC = bsxfun(@times,(1+sign(sM))/2,(f_L + bsxfun(@times,s_minus,(WML - W_L)))) +...
    bsxfun(@times,(1-sign(sM))/2,(f_R + bsxfun(@times,s_plus,(WMR - W_R))));

end

