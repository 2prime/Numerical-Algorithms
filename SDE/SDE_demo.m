clear;
clc;

%%
%Simulator of SDE:
%         dX(t)=aX(t)+bX(t)dW_t,X(0)=1
%
%Author: 2prime
%
%Include:
%     - Euler-Maruyama Secheme
%     - Milstein Secheme
%
%%

%%
%parameter selection

n = 100; a = 5; b=2;

dt = 1/n;
k = a - (b^2)/2;
EM_answer = ones(1,n+1);
Milstein = ones(1,n+1);
true_answer=ones(1,n+1);
plot_base = ones(1,n+1);

white_noise = normrnd(0,sqrt(dt),1,n);


for i = 1:n
    EM_answer(i+1)=EM_answer(i)*(1+a*dt)+b*EM_answer(i)*white_noise(i);
    Milstein(i+1)=Milstein(i)*(1+a*dt)+b*Milstein(i)*white_noise(i)+b*Milstein(i)*(b/2)*(white_noise(i)^2-dt);
    true_answer(i+1)=exp(k*(i/n)+b*sum(white_noise(1:i)));
    plot_base(i+1)=exp(k*(i/n));
end

plot(EM_answer)
hold on
plot(Milstein)
hold on
plot(true_answer)
hold on
plot(plot_base)
hold on
legend('Euler-Maruyana','Milstein','True Answer','With out White Noise')
title('Simulate dX(t)=aX(t)+bX(t)dW_t,X(0)=1')