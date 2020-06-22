clear all;
close all;

m1=1;
m2=1;
I1=0.083;
I2=0.33;
lc1=0.5;
lc2=1;
l1=1;
l2=2;
g=9.81;
q1_dot=0;
q2_dot=0;
q1=-1*pi/2;
q2=0;
q1_dot_ref=0;
q1_ref=pi/2;
kp=16;
kd=8;
angle1=[];
angle2=[];
velocity2=[];
time=[];
t=0;

axis([-6,5,-6,5]);
hold on ;

X2=l1*cos(q1)+l2*cos(q2+q1);
Y2=l1*sin(q1)+l2*sin(q2+q1);
a=[X2,Y2];


A = [    0         0    1.0000         0;
         0         0         0    1.0000;
   12.6292  -12.6926         0         0;
  -14.7488   29.6119         0         0];

B = [    0;
         0;     %calculated in some other files and copied the values which are in accordance with the paper
   -3.0147;
    6.0332];

Q = [1,0,0,0;
     0,1,0,0;
     0,0,1,0;
     0,0,0,1];

R=1;



K =lqr(A,B,Q,R);

angle1 = zeros(600,1);
angle2 = zeros(600,1);
time = zeros(600,1);
velocity2 = zeros(600,1);

%q1_doubledot=v1
%q2_doubledot=v2
tic
for i=1:5000
    
        t=t+0.01;
        d11 =  m1*(lc1^2) + m2*((l1^2)+(lc2^2)+2*l1*lc2*cos(q2))+I1+I2;
        h1  =  -1*m2*l1*lc2*sin(q2)*((q2_dot)^(2)+2*q1_dot*q2_dot); 
        phi1= (m1*lc1+m2*l1)*g*cos(q1) + m2*lc2*g*cos(q1+q2);
        d12 =  m2*(lc2^(2)+l1*lc2*cos(q2))+I2;
%       d1_bar = d21-(d22*d11/d12);
%       h1_bar = h2-(d22*h1/d12);
%       phi1_bar = phi2-(d22*phi1/d12); 
%     
    
        angle1(i,1)= q1;
        angle2(i,1)= q2;
        time(i,1)  = t;
        velocity2(i,1) = q2_dot;
       
        v1=kd*(q1_dot_ref-q1_dot)+kp*(q1_ref-q1);
        v2=(-d11*v1-phi1-h1)/d12;
%         torque = d1_bar*v1+h1_bar+ phi_bar;
        [q1,q1_dot]=finitediff(v1,q1_dot,q1);
        [q2,q2_dot]=finitediff(v2,q2_dot,q2);
        q2= wrapTo2Pi(q2);  
        x=[q1;q2;q1_dot;q2_dot];
        
         
%     end
      
    
    X1=l1*cos(q1);
    Y1=l1*sin(q1);
    X2=l1*cos(q1)+l2*cos(q2+q1);
    Y2=l1*sin(q1)+l2*sin(q2+q1);
    p1= plot([0,X1],[0,Y1],'k');
    hold on;
    p2=plot([X1,X2],[Y1,Y2],'k');
     s1=scatter(X1,Y1,100,'filled','red');
    s2=scatter(X2,Y2,100,'filled','red');
    plot([a(1,1),X2],[a(1,2),Y2],'k');
    a=[X2,Y2];
    pause(0.001);
     delete(s1);
     delete(s2);
    delete(p1);
    delete(p2);
        if (q1>= 1.50 && q2<=0.2 && q2_dot<=0.15)
            display(i);
%           toc;
            break
        end
end

%%
for i=1:600
         display(i);
%         fprintf('hello');
        t=t+0.001;
        x=(x-[pi/2;0;0;0]);
        u=-1*K*x
        x_dot=A*x+B*u;
        
        [q1,q1_dot]= finitediff(x_dot(3,1),x_dot(1,1),q1);
        [q2,q2_dot]= finitediff(x_dot(4,1),x_dot(2,1),q2);
        
        angle1(i,1)=q1;
        angle2(i,1)=q2;
        time(i,1)=t;
        x=[q1;q2;q1_dot;q2_dot];
        
        
        
        
        X1=l1*cos(q1);
        Y1=l1*sin(q1);
        X2=l1*cos(q1)+l2*cos(q2+q1);
        Y2=l1*sin(q1)+l2*sin(q2+q1);
        p1= plot([0,X1],[0,Y1],'k');
        hold on;
        p2=plot([X1,X2],[Y1,Y2],'k');
        s1=scatter(X1,Y1,100,'filled','red');
        s2=scatter(X2,Y2,100,'filled','red');
        plot([a(1,1),X2],[a(1,2),Y2],'k');
        a=[X2,Y2];
        pause(0.01);
        delete(s1);
        delete(s2);
        delete(p1);
        delete(p2);
        
end

toc
close all;
%%
figure(2)
xlabel('time');
ylabel('theta2');
plot(time,angle2);
hold on;
figure(1)
xlabel('time');
ylabel('theta1');
plot(time,angle1);
hold on ;