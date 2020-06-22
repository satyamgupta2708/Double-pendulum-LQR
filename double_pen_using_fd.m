clc
clear all
dt=0.001;
% physical parameters   
m1=1;
m2=1;
l1=1;
l2=2;
lc1=l1/2;
lc2=l2/2;
g=9.81;
tou1=0;
tou2=0;
% main matrix
x_main=[0 0 0 0];
i=1;

%     x=[pi;0.1;9e-14;2e-4];
% 15010
 x=[0; 0; 0; 0];
flag=1;
% loop to calculate x_main matrix
for t=0:dt:30
    
    
    theta1=wrapToPi(x(1));
    theta2=wrapToPi(x(2));
    xdot(1)=x(3);
    xdot(2)=x(4);
 
qdot=[xdot(1); xdot(2)];

% defining damping for pendulum and inertias

I1=0.083;
I2=0.33;
f=[0 0;0 0];

%computing h,c,g 

d=[(m1*(lc1^2)+m2*((l1^2)+2*l1*lc2*cos(theta2)+(lc2^2)))+I1+I2, m2*(lc2^2)+(m2*l1*lc2*cos(theta2))+I2;m2*((lc2^2))+(m2*l1*lc2*cos(theta2))+I2,(m2*(lc2^2))+I2];

hash = -m2*l1*lc2*sin(theta2);
 h = [0,hash*xdot(1)+ hash*(xdot(2)+xdot(1));...
     -hash*xdot(1), 0];
% c=[ 0,-(m2*l1*lc2*sin(theta2)*((xdot(2))))-(2*m2*l1*lc2*sin(theta2)*xdot(1)); (m2*l1*lc2*xdot(1)*sin(theta2)) 0];

grav=[((m1*lc1)+m2*l1)*g*sin(theta1)+(m2*lc2*g*sin(theta1+theta2)),g*m2*lc2*sin(theta1+theta2)]';
% defining parameters for applying pfl
d1_prime=d(2,1)-((d(2,2)*d(1,1))/d(1,2));

h_prime=(h(2,1)*qdot(1))-((d(2,2)*(h(1,1)+h(1,2))*qdot(2))/d(1,2));

g_prime=grav(2,1)-((d(2,2)*grav(1,1))/d(1,2));

kd=8;
kp=16;

v1=0+kd*(0-xdot(1))+(kp*(pi-theta1));
%    if((x(1)~=pi)&&(x(2)~=0.1066)&&(x(3)~=8.9702e-14)&&(x(4)~=-0.0856)&&(flag==1))
% checking conditions for switching between pid and lqr
    if(flag==1) 
    tou2=((d1_prime*v1)+h_prime+g_prime);   %torque for pid
    end

    if((x_main(i,1)<=pi)&&(x_main(i,1)>=3.13)&&((x_main(i,2))<=0.1)&&(x_main(i,2)>=0.0980)&&(x_main(i,3)<=8.999e-14)&&(x_main(i,3)>=0)&&(x_main(i,4)>=-9e-4)&&(x_main(i,4)<=4e-4))
        flag=0;
        disp(10);
    end   
    
%      torque for lqr
     y=x_main(i,:);
     y=y';
    if(flag==0)
          tou2=-([-246.6219,-98.7465,-106.4922,-50.1514]*(y-[pi;0;0;0]));  
    end    
% % choosing torques for applying lqr

%   tou2=-([-120.8777,-49.2658,-53.1548,-21.4757]*(x-[pi;0;0;0])); 
%   
%  else
% parameters for using pd

 
%  end
% torque matrix
tou=[tou1;tou2];
% method to calculate acceleration
qddot=d\(tou-(h*qdot)-grav-(f*qdot));
% assigning variables to 
xdot(3,:)=qddot(1,1);
xdot(4,:)=qddot(2,1);

x(3:4,1)=x(3:4,1)+(xdot(3:4,1)*dt);
x(1:2,1)=x(1:2,1)+(x(3:4,1)*dt);
i=i+1;
x_main(i,1)=wrapToPi(x(1,1));
x_main(i,2)=wrapToPi(x(2,1));
x_main(i,3)=x(3,1);
x_main(i,4)=x(4,1);

% 
end

o=x_main;
axis([-4 4 -4 4] )
for j=1:round(30001*0.1/30):30001
tic;

% links movement
% first link
crank=line([0 sin(o(j,1))],[0 -cos(o(j,1))],'linewidth',2,'MarkerSize',1);
bob1=viscircles([sin(o(j,1)) -cos(o(j,1))],0.05);

% second link
crank2=line([sin(o(j,1)) (sin(o(j,1))+2*sin(o(j,1)+(o(j,2))))],[-cos(o(j,1)) (-cos(o(j,1))-2*cos(o(j,1)+o(j,2)))],'linewidth',2);
 bob2=viscircles([sin(o(j,1))+2*(sin(o(j,1)+o(j,2))) -cos(o(j,1))-2*cos(o(j,1)+o(j,2))],0.05);
pause(0.1-toc)
delete(crank);
delete(crank2);
delete(bob1);
delete(bob2);

toc;

end

