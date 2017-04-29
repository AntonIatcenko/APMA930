%% Pressure Poisson on an Annulus   
% Domain - annulus [1, 2] x [0, 2\pi]
% BC: inhomogeneous Dirichlet in r, periodic in \theta
%% Parameters                       
r1 = 1; r2 = 2;         % domain boundaries
grids = 2.^(2:6);       % trial grids
M = length(grids);      
errors = zeros(1, M);   % preallocating for errors
g = [0 1 0; 0 1 0];     % argument to get the differentiation matrix, see [2] for details
%% Trials                           
parfor j=1:M   % code will work exactly the same with standart for loop, but slower
    %-Grids----------------------------------------------------------------------------------------------------
    Nr = grids(j);  Nth = grids(j); 
    % Radial direction (Neumann BC)
    [rN, pDrr, pDr] = cheb2bc(Nr, g);   rN = -rN*(r2-r1)/2 + (r1+r2)/2;
    % Angular direction (periodic BC)
    [th, Dth] = fourdif(Nth, 1);        [~, Dthth] = fourdif(Nth, 2);
    % Two dimensional grid
    [RR, ThTh] = meshgrid(rN, th);
    %-Operators----------------------------------------------------------------------------------------------- 
    % Base & radial factors
    Irp = eye(Nr);  Ith = eye(Nth); poR = diag(1./rN);  poR2 = poR.^2;
    % Rescaling the operators and reversing the convention
    pDr = -pDr*2/(r2 - r1);    pDrr = pDrr*4/((r2 - r1)^2);
    % Laplacian
    pLap = kron(pDrr, Ith) + kron(poR*pDr, Ith) + kron(poR2, Dthth);
    % Tweaking the pressure Laplacian
    pLap(Nr*Nth/2, :) = 0;  pLap(Nr*Nth/2, Nr*Nth/2) = 1;
    %-Functions------------------------------------------------------------------------------------------------
    u = exp(sin(3*ThTh)).*(RR.^3/3 - 3*RR.^2/2 + 2*RR) + cos(5*RR*pi);
    f = (1./RR).*( exp(sin((3*ThTh))).*(RR.^2 - 3*RR + 2) - 5*pi*sin(5*pi*RR)...
        + RR.*(exp(sin((3*ThTh))).*(2*RR - 3) - 25*pi^2*cos(5*pi*RR))...
        - (9./RR).*sin((3*ThTh)).*exp(sin((3*ThTh))).*(RR.^3/3 - 3*RR.^2/2 ...
        + 2*RR) + (9./RR).*(cos(3*ThTh)).^2.*exp(sin((3*ThTh))).*((RR.^ 3)/3 ...
        - 3*(RR.^2)/2 + (2*RR)));
    fvec = f(:);    uvec = u(:);   fvec(Nr*Nth/2) = uvec(Nr*Nth/2);
    %-Solving--------------------------------------------------------------------------------------------------
    unum = pLap\fvec;  errors(j) = norm(unum-uvec, inf);
    %-Information----------------------------------------------------------------------------------------------
    %fprintf('Grid size = %1.1f x %1.1f, error = %1.5g \n', Nr, Nth, errors(j))
end
%% Convergence results              
figure(1), loglog(grids, errors, '.', 'markersize', 30)
title('Error in Supremum Norm'), xlabel('log(N)'), ylabel('log(error)')