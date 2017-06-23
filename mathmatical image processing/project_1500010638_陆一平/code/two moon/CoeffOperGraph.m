function gamma=CoeffOperGraph(op,alpha,beta,a,b)
[nD Level]=size(alpha);

if strcmp(op,'norm1') || strcmp(op,'norm2')
    gamma=0;
end
for j=1:Level
    for k=1:nD
        if op=='-'
            gamma{k,j}=alpha{k,j}-beta{k,j};
        elseif op=='+'
            gamma{k,j}=alpha{k,j}+beta{k,j};
        elseif strcmp(op,'*c')
            gamma{k,j}=alpha{k,j}*beta;
        elseif strcmp(op,'*+')
            gamma{k,j}=alpha{k,j}*a+beta{k,j}*b;
        elseif op=='h'
            if k==1
                gamma{k,j}=alpha{k,j};
            else
                gamma{k,j}=alpha{k,j}.*(abs(alpha{k,j})>=beta{k,j});
            end
        elseif op=='s'
            if k==1
                gamma{k,j}=alpha{k,j};
            else
                gamma{k,j}=(alpha{k,j}-beta{k,j}.*sign(alpha{k,j})).*(abs(alpha{k,j})>=beta{k,j});
            end
        elseif strcmp(op,'s_band')
            for l=1:length(a)
                if k==a(l)
                    gamma{k,j}=(alpha{k,j}-beta{k,j}.*sign(alpha{k,j})).*(abs(alpha{k,j})>=beta{k,j});
                else
                    gamma{k,j}=alpha{k,j};
                end
            end
        elseif strcmp(op,'norm1')
            gamma=gamma+norm(alpha{k,j},1);
        elseif strcmp(op,'norm2')
            gamma=gamma+norm(alpha{k,j});
        end
    end
end