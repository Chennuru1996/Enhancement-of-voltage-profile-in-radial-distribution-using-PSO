clc;
clear all;
%c=a*P^2+b*P+C(Second degree polynomial); 
%unit no a b c Pmin Pmax
A=[1 0.001562 7.92 561 100 600
2 0.004820 7.97 078 050 200
3 0.001940 7.85 310 100 400]; 
B=[0.0002940 0.0000901 -0.0000507
0.0000901 0.0005210 0.0000953
-0.0000507 0.0000953 0.0000676]/100;
B0=[0.01890 -0.00342 -0.007660]*0.01; 
B00=0.40357;
%load unit3data.m
%A=unit3data; 
pop=50;
npar=3;
XT=585.33; maxit=150;
c1=2.0;
c2=2.0;
wmax=0.9; 
wmin=0.4;     %c2=4-c1;
c3=1; 
 R=10;             %R=10;      %original value 
tol=0.00001;
%----reading input data-----
for i=1:npar
xmin(i)=A(i,5); 
xmax(i)=A(i,6);
a(i)=A(i,2); 
b(i)=A(i,3); 
c(i)=A(i,4);
end
%_______maximum and minimum values for velocities_____ 
for i=1:npar
vmax(i)=(xmax(i)-xmin(i))/R;
vmin(i)=-vmax(i);      %(Xmax(i)-Xmin(i))/R; 
end
%-----random generation of velocities------ 
for pp=1:pop
for ii=1:npar 
 vel(pp,ii)=vmin(ii)+(vmax(ii)-vmin(ii))*rand;
end
end
%------random generation of x i.e.,Pg----- 
for k=1:pop
sumx=0; 
for j=1:npar
x(k,j)=xmin(j)+(xmax(j)-xmin(j))*rand; 
end

%adjustment of losses and loan 
ploss=0;
%sumx=0;
sumx=0;
for j=1:npar 
    X(j)=x(k,j); 
    sumx=sumx+X(j);
end

diff=sumx-XT-ploss; 
while abs(diff)>=tol
for j=1:npar 
    X(j)=X(j)-(diff/npar); 
    if X(j)<xmin(j) 
        X(j)=xmin(j);
elseif X(j)>xmax(j) 
    X(j)=xmax(j);
end
end

L=0;
for j=1:npar
L=L+B0(j)*X(j); end
pl=0;
for j=1:npar
for l=1:npar 
    pl=pl+X(l)*B(j,l)*X(l);
end 
end
ploss=pl+L+B00; %Ploss=ploss+X*B*X'; 
diff=sum(X)-XT-ploss; 
end
%---clculation of fitness and cost functions--- 
sum1=0;
for i=1:npar
%xg(k,i)=X(1,i);%X(i);
sum1=sum1+a(i)*X(1,i)^2+b(i)*X(1,i)+c(i)-abs(xmin(i)-X(1,i)); end

xg(k,:)=X;
err=0;
err=sum(X)-XT-ploss;
cost(k)=sum1; fit(k)=1/(1+cost(k)+(abs(err)*10^4)); 
end
%----initillization of fit ,cost,xg----- 
localfit=fit;
localxg=xg;
localcost=cost;
%-----initial best solution----- 
[globalfit,indx]=max(fit); globalfit=fit(indx); globalcost=cost(indx); globalxg=xg(indx,:);
%-----PSO itration starting------- 
itr=0;
while itr<maxit
itr=itr+1; 
w=wmax-(wmax-wmin)/maxit*itr; 
r1=rand(pop,npar); 
r2=rand(pop,npar);
for k=1:pop
for i=1:npar 
    vel(k,i)=w*vel(k,i)+c1*r1(k,i)*(localxg(k,i)-xg(k,i))+c2*r2(k,i)*(globalxg(1,i)-xg(k,i));
end ;
end;
for k=1:pop for i=1:npar
if vel(k,i)<vmin(i)

vel(k,i)=vmin(i);
elseif vel(k,i)>vmax(i);
vel(k,i)=vmax(i); end
end 
end
xg=xg+vel; for k=1:pop
for i=1:npar
if xg(k,i)>xmax(i)
xg(k,i)=xmax(i); elseif xg(k,i)<xmin(i);
xg(k,i)=xmin(i); end
end ;
end;
for k=1:pop sumx=0;
% for j=1:npar
% x(1,j)=xg(k,j);
%% x(1,j)=xmin(j)+(xmax(j)-xmin(j))*rand; %% sumx=sumx+x(1,j);
%end
%adjustment of losses and load
ploss=0;
%sumx=0;
sumx=0;
for j=1:npar
X(j)=xg(k,j);
sumx=sumx+X(j); end

diff=sumx-XT-ploss; while abs(diff)>=tol
for j=1:npar 
    X(j)=X(j)-(diff/npar); 
    if X(j)<xmin(j) 
        X(j)=xmin(j);
elseif X(j)>xmax(j) 
    X(j)=xmax(j);
end 
end
L=0;
for j=1:npar L=L+B0(j)*X(j);
end
pl=0;
for j=1:npar
for l=1:npar 
    pl=pl+X(j)*B(j,l)*X(l);
end
end
ploss=pl+L+B00; 
diff=sum(X)-XT-ploss; end
xg(k,:)=X;
end
for k=1:pop
cost(k)=0;
for i=1:npar 
    cost(k)=cost(k)+a(i)*xg(k,i)^2+b(i)*xg(k,i)+c(i)-abs(xmin(i)-xg(k,i)); 
end
X=xg(k,:);

err1=0;
L=0;
for j=1:npar
L=L+B0(j)*X(j);
end
pl=0;
for j=1:npar
for l=1:npar
pl=pl+X(j)*B(j,l)*X(l);
end
end
ploss=0;
ploss=pl+L+B00; 
diff=sum(X)-XT-ploss; 
err1=sum(X)-XT-ploss; 
fit(k)=1/(1+cost(k)+(abs(err1)*10^4)); 
end
for k=1:pop
if localfit(k)<fit(k)
localfit(k)=fit(k);
localxg(k,:)=xg(k,:); 
localcost(k)=cost(k);
end
%end
[temp,yy]=max(localfit);
if globalfit<temp
globalfit=localfit(yy); 
globalcost=localcost(yy); 
globalxg=localxg(yy,:); 
end

bcost(itr)=globalcost; 
bfit(itr)=globalfit;
[itr globalxg globalcost globalfit]; 
end
end
subplot(2,1,1)
plot(bcost,'r')
xlabel('number of iterations')
ylabel('local cost')
subplot(2,1,2)
plot(bfit,'b')
xlabel('number of iterations') 
ylabel('best fit')
globalcost
globalxg'
L=0;
for j=1:npar
L=L+B0(j)*globalxg(j); 
end
pl=0;
for j=1:npar
for l=1:npar 
    pl=pl+globalxg(j)*B(j,l)*globalxg(l);
end 
end
ploss=0; 
ploss=pl+L+B00 
sum(globalxg)-XT-ploss
