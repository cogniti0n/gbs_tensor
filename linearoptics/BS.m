function prob = BS(HM, nF, N, flagscreenoutput)
%%% Function to calculate the probability with input Fock state nF
%%% INPUTS
% HM : (1,L) cell, the MPS form of the propagated Fock state
% nF : (1,L) array, the input Fock state
% Dmax : int, the maximum bond dimension
% N : N+1 is the truncation for bosonic Fock space -> corresponds to the photon number!
%%% OUTPUTS
% prob : float, the probability with input nF

% TODO: sanity check
if nargin < 4
    flagscreenoutput = false;
end

L = numel(nF);
tStart = tic;
% 1) construct a Fock state MPS
Mg = cell(1,L);
for itL = (1:L)
    Mg{itL} = zeros(1,1,N+1);
    Mg{itL}(nF(itL)+1) = 1;
end

% 2) contract with HM
H = 1;
for itL = (1:L)
    H = updateLeft(H, 2, Mg{itL}, [], [], HM{itL});
end
prob = abs(H)^2;
tEnd = toc(tStart);
if flagscreenoutput
    fprintf('Boson sampling probability calculated | total time = %6.3f \n', tEnd)
end
end