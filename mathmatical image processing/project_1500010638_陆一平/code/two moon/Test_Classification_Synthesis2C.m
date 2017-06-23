clear
clc

rand('seed',3000);
randn('seed',3000);

load GraphDataSyn_2Circles.mat

N=length(FD)/2;
I1=randperm(N);I0=randperm(N);
Srate=100; % 10% known labels.
Iset1=I1(1:Srate);
Iset0=N+I0(1:Srate);
FD1=FD(Iset1);FD0=FD(Iset0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Extract Decomposition/Reconstruction Masks
FrameType='Haar';
[DFilters RFilters]=ExtractMasks(FrameType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=2; % Dilation scale
n=8; % n-1 = Degree of Chebyshev Polynomial Approximation
Lev=1; % Level of transform
%------------------
% Tuning parameters
%------------------
lambda=0.01; % As in:   lambda*||Wu||_1+1/2||u-f||^2
mu=30*lambda; % Parameter from ADMM. 
%------------------
tol=1e-20; % Tolerance for ADMM. Disabled by default.
maxit=500; % Maximum iterations for ADMM

d=d/max(d); % Normalize the degree vector to [0,1]

M=length(L);
J=log(lambda_max/pi)/log(s)+Lev-1; % Dilation level to start the decomposition
W = @(FD)(GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev));
WT = @(d)(GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev));

[u, nstep]=SplitBregGraphClass2(FD1,FD0,Iset1,Iset0,U1,d,mu,lambda,tol,W,WT,maxit,FD);

FDr=real(u);FDr=double(FDr>=0.5);FDr(Iset1)=1;FDr(Iset0)=0;
Error=sum(abs(FDr-FD))/(M-length(Iset1)-length(Iset0))*100;

subplot(121);scatter(D(1,:),D(2,:),5,FD);title('Ground Truth');axis square;
subplot(122);scatter(D(1,:),D(2,:),5,FDr);title(['Estimated with error = ' num2str(Error) '%.']);axis square;
