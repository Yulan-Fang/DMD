clear all, close all, clc
filename='all-states-history.csv';
All = csvread(filename,1,0,[1,3,10472,42]);%load from 2,4
B=All;
B(:,[1])=[];
C=B;
C([1],:)=[];
s=size(C,1);
d=s/56;
t_days=1:1:d;
D_AK=zeros(1,d);
%% Data X
for i=1:d%20 days
    for j=1:56
        D(j,i)=C(j+(i-1)*56,1);
    end
end
dt=1;
X = D(:,1:end-1);
X2 = D(:,2:end);
[U,S,V] = svd(X,'econ');
%%  Compute DMD (Phi are eigenvectors)
r = 5;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;
lambda = diag(eigs);
omega = log(lambda)/dt;
b = Phi\X(:,1);
mm1 = size(X, 2);
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; 
for iter = 1:mm1 ,
time_dynamics (:,iter )=(b.*exp(omega*t(iter )));
end;
Xdmd = Phi * time_dynamics ;

%%  Plot DMD spectrum
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%% plot DMD constructed data
figure
plot(t,Xdmd(1,end:-1:1))
title('Death Number of AK Constructed by DMD')

for count=1:s %1,1+56,1+112,...
    if (mod(count,56)==1)
    D_AK(1,((d+1)-((count-1)/56+1)))=C(count,1);     
    end
end
hold on
plot(t_days,D_AK)
title('Death Number of AK')
%AK Death

%% prediction
mm2 = size(X, 2);
time_dynamics = zeros(r, mm2);
t = (0:mm2 -1)*dt; 
for iter = 1:mm2 ,
time_dynamics (:,iter )=(b.*exp(omega*t(iter )));
end;
Xdmd2 = Phi * time_dynamics ;
hold on
plot(t,Xdmd2(1,end:-1:1))
title('Death Number of AK Constructed by DMD')