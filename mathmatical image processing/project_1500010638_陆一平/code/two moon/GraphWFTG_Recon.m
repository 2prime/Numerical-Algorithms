function FD_rec=GraphWFTG_Recon(d,L,RFilters,n,s,J,Lev)
r=length(RFilters);
a=pi/2; % Consider the domain of masks as [0, pi]
for j=1:r
    [c_rec{j} ~]=ChebyshevApprox(RFilters{j},n);
end
FD_recl=0;
for l=Lev:-1:1
    for j=1:r
        if l==Lev || j>1
            T0F=d{j,l};
        else
            T0F=FD_rec;
        end
        T1F=(s^(-J+l-1)/a*L)*T0F-T0F;
        djl=1/2*c_rec{j}(1)*T0F+c_rec{j}(2)*T1F;
        for k=3:n
            TkF=(2/a*s^(-J+l-1)*L)*T1F-2*T1F-T0F;
            T0F=T1F;
            T1F=TkF;
            djl=djl+c_rec{j}(k)*TkF;
        end
        FD_recl=FD_recl+djl;
    end
    FD_rec=FD_recl;
    FD_recl=0;
end