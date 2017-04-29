function Er = Navier_Stokes_Annulus_function(Re, Tfinal, r1, r2, b1, b2, Nr, Nth, K)
% Same as Navier_Stokes_Annulus, but returning errors and without plotting.
%% Grids                                  
dt = Tfinal/K;   h = (b2-b1)/(r2-r1);
% Radial direction (Neumann BC)
g = [0 1 0; 0 1 0];   [r, pDrr, pDr] = cheb2bc(Nr, g);  r = -r*(r2-r1)/2 + (r1+r2)/2;
% Angular direction (periodic BC)
[th, Dth] = fourdif(Nth, 1);   [~, Dthth] = fourdif(Nth, 2);
% Two dimensional
[RR, ~] = meshgrid(r(2:end-1), th);           % interior grid
[RRfull, ~] = meshgrid(r, th);            % full grid
% Radial factors
oRvec = 1./RR(:);  oR2vec = oRvec.^2;  R = RRfull(:);
%% Operators                              
% Radial differentiation
[~, Dr] = chebdif(Nr, 2); Dr_in = Dr(2:end-1, 2:end-1, 1); Drr_in = Dr(2:end-1, 2:end-1, 2); Dr = Dr(:,:,1);
% Base & radial factors as operators
Irp = eye(Nr);  Ir = eye(Nr-2);  Ith = eye(Nth);  bigI = kron(Ir, Ith);
poR = diag(1./r); poR2 = poR.^2;  oR = diag(1./r(2:end-1));   oR2 = oR.^2;  bigoR2 = kron(oR2, Ith);
% Rescaling the operators and reversing the convention
Dr_in = -Dr_in*2/(r2 - r1); Drr_in = Drr_in*4/((r2 - r1)^2);  Dr = -Dr*2/(r2 - r1);
pDr = -pDr*2/(r2 - r1);  pDrr = pDrr*4/((r2 - r1)^2);
% Grid operators
DR_in = kron(Dr_in, Ith); oRDR = kron(oR*Dr_in, Ith); DTh = kron(Ir, Dth);      
poRDR = kron(poR*pDr, Ith);  pDTh = kron(Irp, Dth);  DR = kron(Dr, Ith); 
Lap = kron(Drr_in, Ith) + oRDR + kron(oR2, Dthth) - bigoR2;   % Laplacian for diffusion
pLap = kron(pDrr, Ith) + poRDR + kron(poR2, Dthth);           % Laplacian for pressure
pLap(Nr*Nth/2, :) = 0;  pLap(Nr*Nth/2, Nr*Nth/2) = 1;         % tweaking the pressure Laplacian
% Prefactoring
[L, U, Q] = lu(bigI - dt*Lap/Re);        % implicit diffusion
[pL, pU, pQ] = lu(pLap);                 % pressure
%% Steady state solution                  
A = (b2*r2 - b1*r1)/(r2^2 - r1^2); B = (b1/r1 - b2/r2)*r1^2*r2^2/(r2^2 - r1^2);  SSin = A*RR + B./RR;
%% Boundary corrector                     
corr = h*(RR-r1) + b1;  corrIn = corr(:);
%% Initial conditions                     
%bess = @(x) besselj(1, x);  v2d = bess(RR)+SSin;
u2d = 0*RR;  v2d = corr;  u = u2d(:);  v = v2d(:)-corrIn;  % adjusting for corrector 
%% Time integration                       
for t=1:K
    % Explicit advection
    u1 = u - dt*( u.*DR_in*u + oRvec.*(v+corrIn).*DTh*u - oRvec.*(v+corrIn).^2 + 2*oR2vec.*DTh*v/Re );
    v1 = v - dt*( u.*(DR_in*v+h) + oRvec.*(v+corrIn).*DTh*v + oRvec.*u.*(v+corrIn) - 2*oR2vec.*DTh*u/Re );
    % Implicit diffusion
    u2 = U\( L\( Q*u1 ) );  v2 = U\( L\( Q*(v1 + oRvec.*dt*h/Re - dt*oR2vec.*corrIn/Re) ) );  
    % Pressure correction 
    u3 = [0*th; u2; 0*th];  v3 = [0*th; v2; 0*th];  rhs = (u3 + pDTh*v3)./R + DR*u3;
    P = pU\( pL\( pQ*rhs ) );  Pr = DR*P;  Pth = (pDTh*P)./R;
    u = u3 - Pr;  v = v3 - Pth;  u = u(Nth+1:end-Nth);  v = v(Nth+1:end-Nth);
    if sum(isnan([u; v])), warning('Blow up at t = %1.2g', t*dt), break, end   
end
%Er = norm( v + corrIn - exp(-dt*t/Re)*bess(RR(:))-SSin(:), 2)/(Nr*Nth);
%Er = norm( [Pr(1:Nth); Pr(end-Nth+1:end)], 2)/(2*Nth);
Er = norm( [Pth(1:Nth); Pth(end-Nth+1:end)], 2)/(2*Nth);
%Er = norm( [Pr(1:Nth); Pr(end-Nth+1:end)], inf);







