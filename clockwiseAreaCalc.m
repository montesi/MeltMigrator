function A=clockwiseAreaCalc(X,Y)
A=0;
V2=[X(2)-X(1);Y(2)-Y(1)];
for i=3:numel(X);
    V1=V2;
    V2=[X(i)-X(1);Y(i)-Y(1);0];
    C=V1(1)*V2(2)-V1(2)*V2(1);
    A=A+C;
end
A=-A/2;