function F=meltFunctionRJ1981(z,T)

Ga=0.5;
gamma=3; % [degC km^(-1)]
deltaT=600; % [degC]
L=4.2e5; % [J kg^(-1)]
cp=1.05e3; % [J kg^(-1)]
Ts=1100+gamma*z;
epsilon=(1+L/(cp*deltaT))^(-1);
if T<=Ts;
%     Tprime=T-(D-z)*Ga;
    Tprime=T;
else
%     Tprime=Ts+epsilon*(T-(D-z)*Ga-Ts);
    Tprime=Ts+epsilon*(T-Ts);
end
F=(Tprime-Ts)./deltaT;

return