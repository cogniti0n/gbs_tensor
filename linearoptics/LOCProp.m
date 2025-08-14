function [HM, Sent] = LOCProp(nFock, Us, Dmax, opType, flagOutput)
%%% INPUTS
% nFock : an L-dimensional vector indicating the input Fock states
% Us : an L*L unitary matrix indicating the linear optical circuit
% Dmax : int, the maximum bond dimension
% opType : either 'naiive' or 'DMRG'
%%% varargin (optional)
% flagOutput : boolean, set to true if you want to output status (false by default)
%&& OUTPUTS
% HM : an (1,L) cell array representing the propagated MPS
% Sent : the von Neumann entropy at the middle of the MPS

% Sanity check
if (size(Us,1) ~= size(Us,2)) 
    error('ERR: Us needs to be a square matrix!');
end 
if norm(Us*Us'-eye(size(Us,1)),'fro')> 1e-3
    error('ERR: Us needs to be unitary!');
end
if numel(nFock) ~= size(Us,1)
    error('ERR: The shapes of nFock and Us must match!');
end
if (opType~="naiive" && opType~="DMRG")
    error('ERR: The fourth argument "opType" needs to be either "naiive" or "DMRG"');
end

if nargin == 4
    flagOutput = true;
end

if opType == "naiive"
    [HM, Sent] = BS_naiive(nFock, Us, Dmax, flagOutput);
end
if opType == "DMRG"
    [HM, Sent] = BS_DMRG(nFock, Us, Dmax, 4, flagOutput);
end

end

function [HM, Sent] = BS_naiive(nFock, Us, Dmax, flagscreenoutput)
%%% The auxillary function to perform naiive MPO calculations onto the vacuum state
% 0) create relevant variables and local space operators
L = numel(nFock);
N = sum(nFock); % # of photons -> we truncate the bosonic Fock space down to an (N+1)-dimensional space
I = eye(N+1);
A = diag(sqrt((1:N)),1); % annihilation operator

% 1) create the vacuum state MPS
M = cell(1,L);
Mloc = zeros(1,1,N+1);
Mloc(:,:,1) = 1;
M(:) = {Mloc};
M{1} = M{1}(1,:,:);
M{end} = M{end}(:,end,:);

% 2) for each creation operator, contruct the MPO then apply it nFock(m) times
Sent = zeros(L,N);
Mprev = M;
tStart = tic;
for itL1 = (1:L) % current mode (a_m^\dag)

    % construct the MPO for U(a_m^\dag)
    Hs = cell(1,L);
    for itL2 = (1:L)
        Hloc = cell(2,2);
        Hloc(:) = {zeros(size(I))};
        Hloc{1,1} = I;
        Hloc{2,2} = I;
        Hloc{2,1} = conj(Us(itL1,itL2))*A';
        Hs(itL2) = {cell2mat(reshape(Hloc, [1 1 size(Hloc,1), size(Hloc,2)]))};
    end
    Hs{1} = Hs{1}(:,:,end,:);
    Hs{end} = Hs{end}(:,:,:,1);

    % apply the MPO nFock(m) times
    for itN = (1:nFock(itL1))
        Mnow = cell(1,L);
        for itL = (1:L)
            Mnow{itL} = contract(Hs{itL},4,2,Mprev{itL},3,3);

            if itL == 1
                Aleft = 1;
            else
                Aleft = conj(Aright);
            end

            Aright = getIdentity(Mnow{itL},3,Mnow{itL},5);
            % contract isometries
            Mnow{itL} = contract(Aleft,3,[1 2],Mnow{itL},5,[2 4]);            
            Mnow{itL} = contract(Mnow{itL},4,[3 4],Aright,3,[1 2],[1 3 2]);
        end
        Mnow{1} = Mnow{1}(1,:,:);
        Mnow{end} = Mnow{end}(:,end,:);
        
        % truncate bond dimensions by performing a 'round trip' of bringing into canonical forms
        % measure the entanglement entropy in the middle site
        [Mnow, ~] = canonForm(Mnow,L,Dmax,[]);
        [~, S] = canonForm(Mnow,floor(L/2),[],[]);
        Sent(itL1,itN) = -sum(S.^2.*log(S.^2));
        % finally, set the current MPS into the current MPS
        Mprev = Mnow;
    end
end
HM = Mprev;
Sent = reshape(Sent.', 1, []);
Sent = Sent(abs(Sent)>eps);
tEnd = toc(tStart);
if flagscreenoutput
    fprintf('Naiive unitary propagation finished | total time = %6.3f \n', tEnd)
end
end

function [HM, Sent] = BS_DMRG(nFock, Us, Dmax, Nsweep, flagscreenoutput)
% 0) create relevant variables and local space operators
L = numel(nFock);
N = sum(nFock); % # of photons -> we truncate the bosonic Fock space down to an (N+1)-dimensional space
I = eye(N+1);
A = diag(sqrt((1:N)),1); % annihilation operator

% 1) create the vacuum state MPS
M = cell(1,L);
Mloc = zeros(1,1,N+1);
Mloc(:,:,1) = 1;
M(:) = {Mloc};
M{1} = M{1}(1,:,:);
M{end} = M{end}(:,end,:);

% 2) for each creation operator, contruct the MPO then apply it nFock(m) times
Sent = zeros(L,N);
Mprev = M;
tStart = tic;
for itL1 = (1:L) % current mode (a_m^\dag)

    % construct the MPO for U(a_m^\dag)
    Hs = cell(1,L);
    for itL2 = (1:L)
        Hloc = cell(2,2);
        Hloc(:) = {zeros(size(I))};
        Hloc{1,1} = I;
        Hloc{2,2} = I;
        Hloc{2,1} = conj(Us(itL1,itL2))*A';
        Hs(itL2) = {cell2mat(reshape(Hloc, [1 1 size(Hloc,1), size(Hloc,2)]))};
    end
    Hs{1} = Hs{1}(:,:,end,:);
    Hs{end} = Hs{end}(:,:,:,1);
    % apply the MPO nFock(m) times
    for itN = (1:nFock(itL1))
        Mnow = cell(1,L);
        for itL = (1:L)
            Mnow{itL} = contract(Hs{itL},4,2,Mprev{itL},3,3);

            if itL == 1
                Aleft = 1;
            else
                Aleft = conj(Aright);
            end

            Aright = getIdentity(Mnow{itL},3,Mnow{itL},5);
            % contract isometries
            Mnow{itL} = contract(Aleft,3,[1 2],Mnow{itL},5,[2 4]);            
            Mnow{itL} = contract(Mnow{itL},4,[3 4],Aright,3,[1 2],[1 3 2]);
        end
        Mnow{1} = Mnow{1}(1,:,:);
        Mnow{end} = Mnow{end}(:,end,:);
        
        % truncate bond dimensions by performing a 'round trip' of bringing into canonical forms
        % measure the entanglement entropy in the middle site
        [Mnow, ~] = canonForm(Mnow,L,[],[]);
        [~, S] = canonForm(Mnow,floor(L/2),[],[]);
        Sent(itL1,itN) = -sum(S.^2.*log(S.^2));
        % finally, set the current MPS into the current MPS
        Mprev = Mnow;
    end

    Mnew = Mprev;
    for itSweep = (1:Nsweep)
        FrightC = cell(1,L);
        FrightC{end} = 1;
        for it = (L-1:-1:1)
            FrightC{it} = updateLeft(FrightC{it+1}, 2, permute(Mnew{it+1}, [2 1 3]), [], [], permute(Mprev{it+1}, [2 1 3]));
        end
        for it1 = (1:L)
            if it1 > 1
                Fleft = updateLeft(Fleft, 2, Mnew{it1-1}, [], [], Mprev{it1-1});
            else
                Fleft = 1;
            end
            F = contract(Fleft, 2, 2, Mprev{it1}, 3, 1);
            F = contract(F, 3, 2, FrightC{it1}, 2, 2);
            Mnew{it1} = F / sqrt(contract(F, 3, [1 2 3], conj(F), 3, [1 2 3]));
            [Mnew{it1},S,Vd] = svdTr(Mnew{it1}, 3, [1 2], Dmax, []);
            Mnew{it1} = permute(Mnew{it1}, [1 3 2]);
            if it1 < L
                Mnew{it1+1} = contract(diag(S)*Vd, 2, 2, Mnew{it1+1}, 3, 1);
            else
                Mnew{it1} = contract(Mnew{it1}, 3, 2, Vd, 2, 1, [1 3 2]);
            end
        end
    end
    Mprev = Mnew;
end
HM = Mprev;
Sent = reshape(Sent.', 1, []);
Sent = Sent(abs(Sent)>eps);
tEnd = toc(tStart);
if flagscreenoutput
    fprintf('DMRG unitary propagation finished | total time = %6.3f \n', tEnd)
end
end