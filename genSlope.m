function R = genSlope(harmFeat)

a = load(harmFeat);
alog = 20.*log10(a+eps);
[m n]=size(alog);
k=1;
%X1 = [1 2 3];
X2 = [1:30];
slope = zeros(length(X2),1);
for i = 1:m
    y = alog(i,:);
%     Y1 = y(X1);
    Y2 = y(X2);
%     [m1] = linFitt(X1,Y1,length(Y1));
    [m2] = linFitt(X2,Y2,length(Y2));
    slope(i)=m2;
end
sumalog = sum(alog);
% Y3 = sumalog(X1);
Y4 = sumalog(X2);
% [m3] = linFitt(X1,Y3,length(Y3));
[m4] = linFitt(X2,Y4,length(Y4));
aggSlope = m4;

R.SLOPE = slope;
R.AggregateSlope = aggSlope;