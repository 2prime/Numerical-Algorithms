function d=GraphWFTG_Decomp(FD,L,DFilters,n,s,J,Lev)
r=length(DFilters);
for j=1:r
    [c{j} ~]=ChebyshevApprox(DFilters{j},n);
end
a=pi/2; % Consider the domain of masks as [0, pi]
%%%% Fast Tight Frame Decomposition (FTFD)
FD1=FD;
for l=1:Lev
    for j=1:r
        T0F=FD1;
        T1F=(s^(-J+l-1)/a*L)*T0F-T0F;
        d{j,l}=1/2*c{j}(1)*T0F+c{j}(2)*T1F;
        for k=3:n
            TkF=(2/a*s^(-J+l-1)*L)*T1F-2*T1F-T0F;
            T0F=T1F;
            T1F=TkF;
            d{j,l}=d{j,l}+c{j}(k)*TkF;
        end
    end
    FD1=d{1,l};
end