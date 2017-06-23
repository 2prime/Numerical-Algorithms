function [u, nstep]=SplitBregGraphClass2(FD1,FD0,Iset1,Iset0,U1,dd,mu,lambda,tol,W,WT,maxit,forig)
tic;
% u=zeros(size(forig));
u=double(U1<0);
M=length(u);
d=W(u);b=W(u);
normg=CoeffOperGraph('norm2',W(forig));
[r Level]=size(d);
for l=1:Level
    for j=1:r
        Thresh{j,l}=lambda/mu*4^(-l+1)*dd;
    end
end
Iset=[Iset1 Iset0];
Isetc=setdiff(1:M, Iset);
for nstep=1:maxit
    WTdb=WT(CoeffOperGraph('-',d,b));
    u(Isetc)=WTdb(Isetc);
    u(Iset)=1./(dd(Iset)+mu).*([FD1;FD0]+mu*WTdb(Iset));
    u(u<0)=0;u(u>1)=1;
    Wu=W(u);
    d=CoeffOperGraph('s',CoeffOperGraph('+',Wu,b),Thresh);
%     d=CoeffOperGraph('s_band',CoeffOperGraph('+',Wu,b),Thresh,[3:r]);
    deltab=CoeffOperGraph('-',Wu,d);
    b=CoeffOperGraph('+',b,deltab);
    residual=CoeffOperGraph('norm2',deltab)/normg;
    ut=u;ut(Iset1)=1;ut(Iset0)=0;
    error=sum(abs(double(ut>=0.5)-forig))/(M-length(Iset))*100;
    if residual<tol
        break;
    end
    Tm=toc;
    if mod(nstep,maxit)==0
        display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);
    end
end
% display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% display('Program Finished.')
% display(['Step = ' num2str(nstep) '; Residual = ' num2str(residual) '; Error = ' num2str(error) '%; Time Elapsed = ' num2str(Tm)]);