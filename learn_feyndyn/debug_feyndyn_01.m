clear all; close all;

% 1. fundamental constants
kb = 1.3806504e-23;        % Joules / Kelvin
hbar = 1.054571628e-34;    % Joule * seconds

% 2. temperature
temperature = 25;
beta = 1/(kb*temperature);

% 3. system hamiltonian
Omega = hbar*1e12*pi/8;         % Rabi frequency in Joules
H = [0 Omega/2; Omega/2 0];     % System hamiltonian in Joules

% 4. initial denisty matrix
rhoInitial = [0 0; 0 1];

% 5. system coupling matrix
Nbath = 1;                         % Number of baths
systemCouplingMatrix = [0 0; 0 1];   % The system is coupled to the bath via |1X1|
%systemCouplingMatrix=[1 0;0 -1]; % The system is coupled to the bath via sigma_z

% 6. spectral distribution function
alpha = 0.027*pi;               % prefactor. In ps^2 / rad^2
wc = 2.2;                       % cutoff frequency. In rad/ps
dw = 0.01;                      % stepsize for w in J(w). In rad/ps
w = dw:dw:14;                   % must start with dw because coth(0)=infinity. In rad/ps
J = alpha*exp(-(w/wc).^2).*w.^3;% spectral density. In rad/ps

J = J*1e12*hbar;                % spectral distribution function. In Joules

w = w*1e12;                     % w in s-1
dw = dw*1e12;                   % dw in s-1

% 7. numerical parameters for Feynman integral
deltaKmax = 9;                  % number of time steps before memory kernel dies
totalT = 10/1e12;               % total time for the simulation, in seconds
dt = 33/totalT;                 % size of time step (delta t)
finalPoint = round(dt*totalT);  % number of timesteps in total
allPointsORjustFinalPoint = 'allPoints'; % do you just want the density matrix at time=totalT ? or do you want it at every point until then
%cpuORgpu = 'cpu';

% 8. build input array that represents the densitry matrix
% (flattened to a column vector) as a function of time (where the columns represent the time steps)
rho=zeros(numel(H),finalPoint+1);
rho(:,1) = reshape(rhoInitial.',[],1);
% I'm not sure about the transpose, since I'm not sure about
% labeling of the states in the column vector for the purposes of the final sum function.

% 9. run the program and plot only <0|rho|0> and <1|rho|1> as a function of time
wholeDensityMatrixOrJustDiagonals = 'wholeDensityMatrix'; % 'justDiagonals';

%function [rho,elapsedTime] = debug_feyndyn_01(Nbath, ...
%    finalPoint, deltaKmax, totalT, rho, H, systemCouplingMatrix, ...
%    w, dw, J, temperature, wholeDensityMatrixOrJustDiagonals, ...
%    allPointsORjustFinalPoint, cpuORgpu)

% 1. Fundamental constansts
%kb = 1.3806504*10^(-23); % Joules / Kelvin
%hbar = 1.054571628*10^(-34); % Joules * seconds
%beta = 1/(kb*temperature);

% 2. Setup arrays
M = length(H);
M2 = M^2;
Svector = eig(systemCouplingMatrix).';
Svectors = diag(ones(Nbath,1));
diagonals = diag(reshape(1:M2,M,M));
upperTriangle = find(triu(reshape(1:M2,M,M)));
initialPoint = 0;
dt = totalT/finalPoint;

% to determine how big to make the K array.
% This is only necessary for time-dependent OQS hamiltonians.
gridOfTimeIndices = initialPoint:1:finalPoint;

U = expm(-1i*H*dt/hbar);
Uback = U';


% 3. Calculation of eta coefficients

w = [-fliplr(w) w];
J = [-fliplr(J) J];
temp = beta*hbar*w/2;
MakMak = (1/(hbar*pi))*J.*exp(temp)./(sinh(temp).*w.^2);

% the -1 is because eta1 isn't meant to go up to deltaKmax,
% that's what eta 5 is for
for deltaK=1:deltaKmax%-1
    eta1(deltaK) = sum(2*MakMak.*((sin(w*dt/2)).^2).*exp(-1i*w*dt*deltaK))*dw;
end

eta2 = sum(0.5*MakMak.*(1-exp(-1i*w*dt)))*dw;
for N = 2:deltaKmax+1
    % N=1 is a dummy index for the eta_00 case, which is covered by eta4 and eta2
    eta3(N) = sum(2*MakMak.*((sin(w*dt/4)).^2).*exp(-1i*w*((N-1)*dt - dt/2)))*dw;
    % the N-1 is because, N is actually starting at 2 and going to deltaKmax+1,
    % instead of going from 1 to deltaKmax
end

eta4 = sum(0.5*MakMak.*(1-exp(-1i*w*dt/2)))*dw;
% k=1 is a dummy index for the eta_00 case, which is covered by eta4 and eta2
for k=2:deltaKmax+1
    % the k-1 is because, k is actually starting at 2 and going to deltaKmax+1,
    % instead of going from 1 to deltaKmax
    eta5(k) = sum(2*MakMak.*sin(w*dt/4).*sin(w*dt/2).*exp(-1i*w*((k-1)*dt-dt/4)))*dw;
end

%-1 %because eta6 isn't meant to go up to deltaKmax, that's what eta 3 is for
for NminusK = 1:deltaKmax
    eta6(NminusK) = sum(2*MakMak.*sin(w*dt/4).*sin(w*dt/2).*exp(-1i*w*((NminusK)*dt-dt/4)))*dw;
end

etaK = zeros(1,deltaKmax+1);
etaK(1) = eta2;
etaK(2:deltaKmax+1) = eta1; %etaK(deltaKmax+1)=eta5(deltaKmax);

etaN = zeros(1,deltaKmax+1);
etaN(1) = eta4;
etaN(2:deltaKmax+1) = eta6;%etaN(deltaKmax+1)=eta3(deltaKmax);
kernel = etaK;
kernelConj = conj(kernel);
kernelEnd = etaN;
kernelEndConj = conj(kernelEnd);
kernelTD = eta5;
kernelTDConj = conj(kernelTD);
kernelTDend = eta3;
kernelTDendConj = conj(kernelTDend);
%kernelTemp(1)=eta5_10, (2)=eta5_20, etc.. kernelTempEnd is same but with eta3

kernelTD(1) = eta4;
kernelTDend(1) = eta4;


%% 4. Calculation of K tensor and I tensors
% The way that the indices are defined for the influence functional
indicesForF = my_npermutek(1:M,4);
indicesForK = [indicesForF(:,3) , indicesForF(:,1) , indicesForF(:,2) ,indicesForF(:,4)];

%row indicesForK(:,1) of U and column indicesForK(:,2) of U .. see page 8 of lab book
K = U(indicesForK(:,1) + M*(indicesForK(:,2)-1)) .* Uback(indicesForK(:,3) + M*(indicesForK(:,4)-1));
% for time dependent OQS hamiltonians
K = repmat(K,1,length(gridOfTimeIndices));

[ia,ib] = ndgrid(1:M);

% would not need to transpose the things in the next three lines
% if Mvector was already a column vector (1:M)'
Sdiff = (kron(Svector,ones(1,M))-kron(ones(1,M),Svector))';

%eq(14) on pg 4601 of tensor prop. I. repeated values are for n=1, changing values are for n=0
I_00 = my_expand(exp(-Sdiff.*(kernelEnd(1)*Svector(1,ib) - kernelEndConj(1)*Svector(1,ia)).'),[M2,1]);
I_mk = zeros(M2^2,deltaKmax+1);
I_mkEnd = I_mk;
I_mkTD = I_mk;
I_mkTDend = I_mk; %I_00 means I0(sk), I_mk means I(sk,s minus k) ?
%eq(14) on pg 4601 of tensor prop. I. repeated values are for 0 (or just n-1?),
% changing values are for n %might be better to rearrange stuff so that you can
% use kron if repmat is slow in fortran
% % is it better computationally to do this in one line or separate calcs ?
% This is delatK=0, ie, eq 14 on pg 4604 of Tensor Prop. I, the diff btwn this
% and the line above is that this is for sk>=1 while line above is for sk=0 ? eta+kk ?
I_mk(:,1) = kron(ones(M2,1),exp(-Sdiff.*(kernel(1)*Svector(1,ib)-kernelConj(1)*Svector(1,ia)).'));
I_mkEnd(:,1) = kron(ones(M2,1),exp(-Sdiff.*(kernelEnd(1)*Svector(1,ib)-kernelEndConj(1)*Svector(1,ia)).')); %repeated values are for 0, changing values are for n %might be better to rearrange stuff so that you can use kron if repmat is slow in fortran, %is it better computationally to do this in one line or separate calcs ? This is delatK=0, eta_NN ?

%should go up to deltaKmax-1  and deltaKmax should be treated separately
for deltaK = 1:deltaKmax
    I_mk(:,1+deltaK)=exp(-kron((kernel(1+deltaK)*Svector(1,ib)-kernelConj(1+deltaK)*Svector(1,ia)).',Sdiff));
    %can't we just make kernel a vector and use kernel(deltaK) ? , %saving kernelConj reduces number
    % of times 'conj.m' has to be called, but adds another variable, which is worse ?
    % %is it bad to have such a big argument in kron ? or should i save it and make 2 lines,
    % depending on your units for kernel, you might need to factor  the argument by 1/hbar.
    % The reshaping afterwards must be redundant, we must be able to use kron in such a way that
    % the answer ends up in matrix form ,
    I_mkTD(:,1+deltaK)=exp(-kron((kernelTD(1+deltaK)*Svector(1,ib)-kernelTDConj(1+deltaK)*Svector(1,ia)).',Sdiff));
    I_mkEnd(:,1+deltaK)=exp(-kron((kernelEnd(1+deltaK)*Svector(1,ib)-kernelEndConj(1+deltaK)*Svector(1,ia)).',Sdiff));
    I_mkTDend(:,1+deltaK)=exp(-kron((kernelTDend(1+deltaK)*Svector(1,ib)-kernelTDendConj(1+deltaK)*Svector(1,ia)).',Sdiff));
end %kernel(1)=eta_kk, kernel(2)=eta_{k+1,k}

%% 4.1 Calculation of I tensors with N identical baths
for jj = 2:Nbath
    Svector2 = Svectors(Nbath+1-jj,:); %relationship between Svector and Svector2
    Sdiff2 = (kron(Svector2,ones(1,M))-kron(ones(1,M),Svector2))';
    Ij_00 = expand(exp(-Sdiff2.*(kernelEnd(1)*Svector2(1,ib)-kernelEndConj(1)*Svector2(1,ia)).'),[M2,1]);
    Ij_mk = zeros(M2^2,deltaKmax+1);
    Ij_mkEnd = Ij_mk;
    Ij_mkTD = Ij_mk;
    Ij_mkTDend = Ij_mk;
    Ij_mk(:,1) = kron(ones(M2,1),exp(-Sdiff2.*(kernel(1)*Svector2(1,ib)-kernelConj(1)*Svector2(1,ia)).'));
    Ij_mkEnd(:,1) = kron(ones(M2,1),exp(-Sdiff2.*(kernelEnd(1)*Svector2(1,ib)-kernelEndConj(1)*Svector2(1,ia)).'));

    for deltaK=1:deltaKmax
        Ij_mk(:,1+deltaK)=exp(-kron((kernel(1+deltaK)*Svector2(1,ib)-kernelConj(1+deltaK)*Svector2(1,ia)).',Sdiff2));
        Ij_mkTD(:,1+deltaK)=exp(-kron((kernelTD(1+deltaK)*Svector2(1,ib)-kernelTDConj(1+deltaK)*Svector2(1,ia)).',Sdiff2));
        Ij_mkEnd(:,1+deltaK)=exp(-kron((kernelEnd(1+deltaK)*Svector2(1,ib)-kernelEndConj(1+deltaK)*Svector2(1,ia)).',Sdiff2));
        Ij_mkTDend(:,1+deltaK)=exp(-kron((kernelTDend(1+deltaK)*Svector2(1,ib)-kernelTDendConj(1+deltaK)*Svector2(1,ia)).',Sdiff2));
    end %kernel(1)=eta_kk, kernel(2)=eta_{k+1,k}

    I_00 = I_00.*Ij_00;
    I_mk = I_mk.*Ij_mk;
    I_mkTD = I_mkTD.*Ij_mkTD;
    I_mkEnd(:,1) = I_mkEnd(:,1).*Ij_mkEnd(:,1);
    I_mkTDend = I_mkTDend.*Ij_mkTDend;
end

%% 5. Propagation of rho for the first deltaKmax timesteps
A = K(:,1) .* I_mkTD(:,1+1) .* I_mk(:,1) .* I_00 .* my_expand(rho(:,1),[M2,1]);
%I_mkTD=I_k0,k=1  I_mk=I_kk, k=1 , I_00=I_NN,k=0
Aend = K(:,1) .* I_mkTDend(:,1+1) .* I_mkEnd(:,1) .* I_00 .* my_expand(rho(:,1),[M2,1]);
%I_mkTD=I_N0,N=1  I_mkEnd=I_NN, k=1 , I_00=I_NN,N=0
%K(:,1)=[]; % only for time-dependent system Hamiltonians

rho(:,2) = sum(reshape(Aend,M2,M2).');
Aend = [];

indices = uint32(my_npermutek(cast(1:M2,'single'),1+1));
%makes 1:4 single precision. One can make them int(8) but then you'd have to
% use Jan simons's thing (or you could use Matt Fig's MEX, since you're storing it anyway)

%could go to deltaKmax-1 and incorporate the deltaKmax case into the next forloop,
% but the former way uses I_mkTD and the latter does not.
for J = 2:deltaKmax
    % Making 1:M single precision might help
    indices = horzcat(my_expand(indices,[M2,1]),repmat((1:M2)',size(indices,1),1));
    A = (my_expand(A,[M2,1]));
    %A=A.*repmat(K_k,M2^(J-1),1); % or just read indices like in below lines, whatever's faster [this is implemented in line below because it's needed for OFPF]
    A = A.*K((indices(:,end-1)-1)*M2+indices(:,end),1); % %Or just make column of K = J, and forget about the deleting of columns.
    K(:,1) = [];
    Aend = A;
    for k = 0:J-1
        % this shouldn't be redone each time, since these lines have no dependence on J
        A = A.*I_mk((indices(:,end-k)-1)*M2+indices(:,end),k+1);
        Aend = Aend.*I_mkEnd((indices(:,end-k)-1)*M2+indices(:,end),k+1);
    end
    %the k=J case for forloop above. uses eta_k0
    A = A.*I_mkTD((indices(:,end-J)-1)*M2+indices(:,end),J+1);
    % the k=J case for forloop above. uses eta_N0
    Aend = Aend.*I_mkTDend((indices(:,end-J)-1)*M2+indices(:,end),J+1);
    %whos('A','Aend','indices')

    %weakPaths=find(abs(A)<tol); %in 2001 Sim, it's <=
    %A(weakPaths)=0;indices(weakPaths,:)=0; % might be better to just remove these rows completely, and then calculate rho using information in the leftover INDICES array. If you remove rows of A, you have to also remove rows of indices to allow A.*I_mk(indices) to work
    %A=sparse(A);indices=sparse(indices); %now there's a problem with indices=horzcat.... above. Not only does size(indices,1) have to be changed into length(find(indices(:,1)), but expand will also expand the 0's even though you don't want that.
    %A(weakPaths)=[];indices(weakPaths,:)=[];weakPaths=[]; %might be better to just remove these rows completely, and then calculate rho using information in the leftover INDICES array. If you remove rows of A, you have to also remove rows of indices to allow A.*I_mk(indices) to work
    %rho(:,J+1)=sum(reshape(A,M2,M2^(J)).'); %watch out, ' is conjugate transpose, we only want to transpose
    for Q=1:size(rho,1)
        rho(Q,J+1) = sum(Aend(indices(:,end)==Q)); %do as array function ?
    end
end

%% 6. Propagation of rho for the timesteps after deltaKmax
tic;
if finalPoint > deltaKmax

    KtimesIpermanent = K((indices(:,end-1)-1)*M2+indices(:,end),1); %predefine ?
    KtimesIpermanentEnd = KtimesIpermanent;
    %whos
    for k=0:deltaKmax % it's probably faster if this tensor is just built as A is in the above forloop ... in fact it seems reading indices is slow, so perhaps even the previous iterations should be done like the tensor propagator below:
        KtimesIpermanent = KtimesIpermanent.*I_mk((indices(:,end-k)-1)*M2+indices(:,end),k+1); %when k=deltaKmax, we're using I_mk(:,1+deltaKmax) which used kernel(1+deltaK) which uses eta 5
        KtimesIpermanentEnd = KtimesIpermanentEnd.*I_mkEnd((indices(:,end-k)-1)*M2+indices(:,end),k+1);%when k=deltaKmax, we're using I_mkEnd(:,1+deltaKmax) which used kernelEnd(1+deltaK) which uses eta 3
    end
    indices=[]; %clear('indices') %when not using parfor
    KtimesIpermanent = reshape(KtimesIpermanent,M2,[]);
    KtimesIpermanentEnd = reshape(KtimesIpermanentEnd,M2,[]);

    for J = deltaKmax+1:finalPoint
        Aend = reshape(bsxfun(@times,reshape(sum(reshape(A,[],M2),2),1,[]),KtimesIpermanentEnd),[],1); %see profile.txt for discussion about speed
        A = reshape(bsxfun(@times,reshape(sum(reshape(A,[],M2),2),1,[]),KtimesIpermanent),[],1); %could be sped up by using some of the calculations in line above

        Aend = reshape(Aend,M2,[]); % for some reason, switching the M2 and [] and removing the ,2 in the line below alters the beginning
        if or(strcmpi(allPointsORjustFinalPoint,'allPoints'),and(strcmpi(allPointsORjustFinalPoint,'justFinalPoint'),J==finalPoint));
            switch wholeDensityMatrixOrJustDiagonals
                case 'justDiagonals'
                    rho(diagonals(1:M-1),J+1)=sum(Aend(diagonals(1:M-1),:),2); % all diagonals but the last
                    rho(diagonals(end),J+1)=1-sum(rho(:,J+1));                 % last diagonal obtained by trace(rho)=1
                otherwise
                    rho(upperTriangle(1:end-1),J+1)=sum(Aend(upperTriangle(1:end-1),:),2); % all of the upper-right triangle of the matrix, except for the last diagonal
                    rho(diagonals(end),J+1)=1-sum(rho(diagonals(1:end-1),J+1));            % last diagonal obtained by trace(rho)=1
            end
            disp(['Time step ' num2str(J) '/' num2str(finalPoint) ' has completed successfully after ' num2str(round(toc)) ' seconds! L=' num2str(length(find(A)))]);
        end
    end
end

elapsedTime = toc;

% assign output ????

%figure(1);hold('on')
%plot(0:totalT/finalPoint:totalT,real(rho(1,:)));
%plot(0:totalT/finalPoint:totalT,real(rho(4,:)));
%xlabel('time (seconds)');

figure(2);hold('on')
plot(0:totalT/finalPoint:totalT,real(rho(1,:)));
plot(0:totalT/finalPoint:totalT,real(rho(3,:)),'r');
plot(0:totalT/finalPoint:totalT,imag(rho(3,:)),'--r');
plot(0:totalT/finalPoint:totalT,real(rho(3,:)),'r');
plot(0:totalT/finalPoint:totalT,-imag(rho(3,:)),'--r');
plot(0:totalT/finalPoint:totalT,real(rho(4,:)));
xlabel('time (seconds)');
