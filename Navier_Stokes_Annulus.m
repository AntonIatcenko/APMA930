%% Navier-Stokes Equations on an Annulus  
% with Chorin pressure correction
% u_t = - uu_r - vu_th/r + v^2/r + ( (ru_r)_r/r + u_{thth}/(r^2) - u/(r^2) - 2v_th/(r^2))/Re - p_r
% v_t = - uv_r - vv_th/r + uv/r + ( (rv_r)_r/r + v_{thth}/(r^2) - v/(r^2) + 2v_th/(r^2))/Re - p_th/r
% (ru_r)_r + v_th = 0
% Domain - annulus [r1, r2] x [0, 2\pi]
% BC: inhomogeneous Dirichlet in r, periodic in \theta
clear, tic
%% Physical Parameters                    
Re = 1e+1;             % Reynolds number
Tfinal = 1;            % final time
r1 = 1; r2 = 2;        % domain boundaries
%bess = @(x) besselj(1, x);
%r1 = fzero(bess, 3.83); r2 = fzero(bess, 16.47);
b1 = 5; b2 = 0;        % boundary conditions 
h = (b2-b1)/(r2-r1);   % slope of the boundary correction
%% Computational Parameters               
Nr = 32;      Nth = 32;                 % resolution in radial and angular directions              
K = 1e+3;     dt = Tfinal/K;            % number and size of time steps
plotgap = 1e+1; numplots = K/plotgap;   % number of time steps between plots and total number of plots
%% Grids                                  
% Radial direction (Neumann BC)
g = [0 1 0; 0 1 0];   [r, pDrr, pDr] = cheb2bc(Nr, g);  r = -r*(r2-r1)/2 + (r1+r2)/2;
% Angular direction (periodic BC)
[th, Dth] = fourdif(Nth, 1);   [~, Dthth] = fourdif(Nth, 2);
% Two dimensional
[RR, ThTh] = meshgrid(r(2:end-1), th);            % interior grid
[RRfull, ThThfull] = meshgrid(r, th);             % full grid
XX = [r'; RRfull].*cos([ThThfull; th'*0+2*pi]);   % plotting grid
YY = [r'; RRfull].*sin([ThThfull; th'*0+2*pi]);   % plotting grid
% Radial factors
oRvec = 1./RR(:);  oR2vec = oRvec.^2;  R = RRfull(:);
% Preallocating for the solution
Udata=zeros(numplots+1, (Nr-2)*Nth); Vdata=zeros(numplots+1, (Nr-2)*Nth); Pdata=zeros(numplots+1, Nr*Nth);
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
corr = h*(RRfull-r1) + b1;  corrext = [corr; corr(1, :)];  corrIn = h*(RR(:)-r1) + b1;
%% Initial conditions                     
u2d = 0*RR; v2d = 0*ThTh; u = u2d(:); v = v2d(:)-corrIn; % adjusting for corrector 
Udata(1, :) = u;  Vdata(1, :) = v+corrIn;   %v2d=10*bess(RR)+SSin;
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
    if mod(t, plotgap) == 0
        Udata(t/plotgap+1, :) = u;  Vdata(t/plotgap+1, :) = v+corrIn;  Pdata(t/plotgap+1, :) = P/dt;
    end
    if sum(isnan([u; v])), warning('Blow up at t = %1.2g', t*dt), break, end   
end, disp('Done with time integration, creating plots.'), toc
%% Plotting                               
vidObj = VideoWriter('Video_TC.avi');  vidObj.FrameRate = 10;  open(vidObj);
% Adding boudary conditions
UwBC = [zeros(numplots+1, Nth) Udata zeros(numplots+1, Nth)];
VwBC = [ones(numplots+1, Nth)*b1 Vdata ones(numplots+1, Nth)*b2];
energies = zeros(2, numplots+1); 
warning('off', 'MATLAB:contour:ConstantData')
for t = 1:numplots+1
    U = reshape(UwBC(t, :), Nth, Nr);   Uext = [U; U(1, :)];
    V = reshape(VwBC(t, :), Nth, Nr);   Vext = [V; V(1, :)];
    P = reshape(Pdata(t, :), Nth, Nr);  Pext = [P; P(1, :)];
    
    fig=figure(1); 
    sp1 = subplot(2, 2, 1); p = get(sp1, 'position'); set(sp1, 'position', p.*[1 1 1.15 1.15]-[0 .03 0 0]);
    contourf(XX, YY, Uext, 10, 'edgecolor', 'none')
    title(['Radial velocity u at time ', num2str(dt*(t-1)*plotgap)]), colorbar, 
    xlabel('X'), ylabel('Y'), axis square
    
    sp2 = subplot(2, 2, 2); p = get(sp2, 'position'); set(sp2, 'position', p.*[1 1 1.15 1.15]-[0 .03 0 0]);
    contourf(XX, YY, Vext, 20, 'edgecolor', 'none')
    title(['Angular velocity u at time ', num2str(dt*(t-1)*plotgap)]), colorbar, 
    xlabel('X'), ylabel('Y'), axis square
    
    sp3 = subplot(2, 2, 3); p = get(sp3, 'position'); set(sp3, 'position', p.*[1 1 1.15 1.15]-[0 .03 0 0]);
    contourf(XX, YY, Pext, 30, 'edgecolor', 'none')
    title(['Pressure at time ', num2str(dt*(t-1)*plotgap)]), colorbar, 
    xlabel('X'), ylabel('Y'), axis square
    
    sp4 = subplot(2, 2, 4); p = get(sp4, 'position'); set(sp4, 'position', p.*[1 1 1.15 1.15]-[0 .03 0 0]);
    %plot(r, 10*exp(-dt*(t-1)*plotgap/Re)*bess(r)+[b1 SSin(Nth/2, :) b2]', 'linewidth', 1), hold on
    plot(r, [b1 SSin(Nth/2, :) b2]', 'linewidth', 1), hold on
    plot(r, Vext(Nth/2, :), '.', 'markersize', 20), xlim([r1 r2]), ylim([-1 b1+1])
    title(['Slice \theta = \pi at time ', num2str(dt*(t-1)*plotgap)])
    %legend('Analytical Solution', 'Numerical Solution'), hold off
    legend('Steady State', 'Numerical Solution'), hold off
    drawnow, currFrame = getframe(fig);   writeVideo(vidObj, currFrame);
    
    energies(1, t) = norm(U, 2)/(Nr*Nth);  energies(2, t) = norm(V, 2)/(Nr*Nth);
end, warning('on', 'MATLAB:contour:ConstantData'), close(vidObj);
figure(2)
subplot(2, 1, 1), plot((0:numplots)*dt*plotgap, energies(1, :), 'linewidth', 3), title('L_2 norm of u')
subplot(2, 1, 2), plot((0:numplots)*dt*plotgap, energies(2, :), 'linewidth', 3), title('L_2 norm of v')