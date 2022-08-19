%--------------------------------------------------------------------------
%-------------- Double Pendulum--------------------------------------------
%--------------------------------------------------------------------------
clc
close all
clear all


%---------Parameters------------------------------------------------------
l1=1; l2=1 ; m1=1 ; m2=1; g=9.81;
%---------initial condition-----------------------------------------------
tspan=25;
theta1=pi/2;
theta1_prime=0;
theta2=pi/2;
theta2_prime=0;
y0=[theta1 theta1_prime theta2 theta2_prime];
[t,y]=ode45(@pend, [0 ,tspan], y0);
%---position of mass 1 and mass 2----------------------------------------
x1=l1*sin(y(:,1));
y1=-l1*cos(y(:,1));
x2=l1*sin(y(:,1))+l2*sin(y(:,3));
y2=-l1*cos(y(:,1))-l2*cos(y(:,3));
%------visualizing the result---------------------------------------------
   figure(1)
   plot(x1,y1,'linewidth',1.5)
   hold on
   plot(x2,y2,'r','linewidth',1.5)
   h=gca; 
   get(h,'fontSize') 
   set(h,'fontSize',14)
   xlabel('X','fontSize',14);
   ylabel('Y','fontSize',14);
   title('Double Pendulum Mass Motion Trace','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   
theta1=pi/2+0.001;
theta1_prime=0;
theta2=pi/2+0.001;
theta2_prime=0;
y0=[theta1 theta1_prime theta2 theta2_prime];
[t,y]=ode45(@pend, [0 ,tspan], y0);
%---position of mass 1 and mass 2----------------------------------------
x1=l1*sin(y(:,1));
y1=-l1*cos(y(:,1));
x2=l1*sin(y(:,1))+l2*sin(y(:,3));
y2=-l1*cos(y(:,1))-l2*cos(y(:,3));
%------visualizing the result---------------------------------------------
   figure(1)
   plot(x1,y1,'linewidth',1.5)
   hold on
   plot(x2,y2,'r','linewidth',1.5)
   h=gca; 
   get(h,'fontSize') 
   set(h,'fontSize',14)
   xlabel('X','fontSize',14);
   ylabel('Y','fontSize',14);
   title('Double Pendulum Mass Motion Trace','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
%{
   figure(2)
   plot(y(:,1),'linewidth',1.5)
   hold on
   plot(y(:,3),'r','linewidth',1.5)
   h=gca; 
   get(h,'fontSize') 
   set(h,'fontSize',14)
   legend('\theta_1','\theta_2')
   xlabel('time','fontSize',14);
   ylabel('theta','fontSize',14);
   title('\theta_1(t=0)=\pi/2 and \theta_2(t=0)=\pi/2','fontsize',14)
   fh = figure(2);
   set(fh, 'color', 'white'); 
   
   figure(3)
   plot(t,x1,'linewidth',1.5)
   hold on
   plot(t,x2,'r','linewidth',1.5)
   h=gca; 
   get(h,'fontSize') 
   set(h,'fontSize',14)
   legend('x1','x2')
   xlabel('Time','fontSize',14);
   ylabel('X','fontSize',14);
   title('Mass x position vs time','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   
   figure(4)
   plot(t,y1,'linewidth',1.5)
   hold on
   plot(t,y2,'r','linewidth',1.5)
   h=gca; 
   get(h,'fontSize') 
   set(h,'fontSize',14)
   legend('y1','y2')
   xlabel('Time','fontSize',14);
   ylabel('Y','fontSize',14);
   title('Mass y position vs time','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   
   
   %----movie of double pendulum--------------------------------------
   
   figure(3)
   Ncount=0;
   fram=0;
   
     for i=1:length(y)
         Ncount=Ncount+1;
         fram=fram+1;
         plot(0, 0,'.','markersize',20);
         hold on
         plot(x1(i),y1(i),'.','markersize',20);
         plot(x2(i),y2(i),'.','markersize',20);
         hold off
         line([0 x1(i)], [0 y1(i)],'Linewidth',2);
         axis([-(l1+l2) l1+l2 -(l1+l2) l1+l2]);
         line([x1(i) x2(i)], [y1(i) y2(i)],'linewidth',2);
         h=gca; 
         get(h,'fontSize') 
         set(h,'fontSize',12)
         xlabel('X','fontSize',12);
         ylabel('Y','fontSize',12);
         title('Chaotic Motion','fontsize',14)
         fh = figure(3);
         set(fh, 'color', 'white'); 
         F=getframe;
         end
      movie(F,fram,20)
   
%}
      
function [yprime] = pend(t, y)
l1=1; l2=2 ; m1=2 ; m2=1; g=9.8;
y_prime=zeros(4,1);
a = (m1+m2)*l1 ;
b = m2*l2*cos(y(1)-y(3)) ;
c = m2*l1*cos(y(1)-y(3)) ;
d = m2*l2 ;
e = -m2*l2*y(4)* y(4)*sin(y(1)-y(3))-g*(m1+m2)*sin(y(1)) ;
f = m2*l1*y(2)*y(2)*sin(y(1)-y(3))-m2*g*sin(y(3)) ;
yprime(1) = y(2);
yprime(3)= y(4) ;
yprime(2)= (e*d-b*f)/(a*d-c*b) ;
yprime(4)= (a*f-c*e)/(a*d-c*b) ;
yprime=yprime';
end
      
 %-----------------------------------------------------------------------     