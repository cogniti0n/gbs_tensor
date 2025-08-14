function prob = GBS(HM, zSq, N, flagscreenoutput)
% Function to calculate the probability with input squeezed state zSq
%%% INPUTS
% HM : (1,L) cell, the MPS form of the propagated Fock state
% zSq : (1,L) array, the squeezing parameters
% N : N+1 is the truncation for bosonic Fock space -> corresponds to the photon number!
%%% OUTPUTS
% prob : float, the probability with input squeezed state

L = numel(HM); % # of sites
if nargin < 5
    flagscreenoutput = false;
end

tStart = tic;
% 1) construct a Gaussian state MPS
Mg = cell(1,L);
for itL = (1:L)
    Mloc = zeros(1,1,N+1);
    Mloc(1,1,:) = squeezed(zSq(itL),N);
    Mg{itL} = Mloc;
end

% 2) contract the MPS with the HM MPS
H = 1;
for itL = (1:L)
    H = updateLeft(H, 2, Mg{itL}, [], [], HM{itL});
end

% 3) calculate the probability
prob = abs(H)^2;
tEnd = toc(tStart);
if flagscreenoutput
    fprintf('Gaussian boson sampling probability calculated | total time = %6.3f \n', tEnd)
end
end