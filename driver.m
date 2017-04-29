%% Driver for Navier_Stokes_Annulus_function
% Runs the function on with different parameters and produces converence
% information
tic
%% Parameters                    
Re = 1e+1;                   % Reynolds number
Tfinal = 1e+3*2.^(5:20);         % final time
%bess = @(x) besselj(1, x);
%r1 = fzero(bess, 3.83);      % domain boundaries
%r2 = fzero(bess, 16.47);
r1 = 1; r2 = 3;
K = 1;
grids = 6:2:32;
M = length(Tfinal);
errors = zeros(1, M);

parfor j=1:M

              % Navier_Stokes_Annulus_function(Re, Tfinal, r1, r2, b1, b2, Nr, Nth, K)
    errors(j) = Navier_Stokes_Annulus_function(Re, 1/Tfinal(j), r1, r2, 5, 0, 32, 32, K)


end
% %% Convergence results              
% figure(1), loglog(Tfinal, errors(1, :), '.', 'markersize', 20)
% title('Error in 2-Norm with b_1 = 0'), xlabel('log(N)'), ylabel('log(error)')
% figure(2), loglog(Tfinal, errors(2, :), '.', 'markersize', 20)
% title('Error in 2-Norm with b_1 = 10'), xlabel('log(N)'), ylabel('log(error)')

%% Convergence results           
p = polyfit( log(Tfinal), log(errors), 1);

figure(1)
loglog(Tfinal, exp(polyval(p, log(Tfinal))), 'linewidth', 2)
hold on
loglog(Tfinal, errors, '.', 'markersize', 20)
hold off
title('Error in 2-Norm with b_1 = 5')
legend(['Rate = ', num2str(-p(1))])
xlabel('log(1/\Delta t)'), ylabel('log(error)')
toc




