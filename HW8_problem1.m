
% Solution for HWk 8 H.E/ Diffusion equation
% part one only
%will print out a 2-d graph that is 'animated' to show how closed-form 
%and how Crank-Nicholson change w.r.t. and how they compare w/ each other
clear all 
clc


%defining parameters in this space
L=pi;
D=.1;
T=10;
F=0;
time =0;
space=0;

%will divide interval into steps equal  INNER space steps
steps=40;

dx=L/(steps+1);
x=0:dx:pi;

%time interval  steps
dt=T/200;
t=0:dt:T;

% defining 'lambda' to help solve 
lambda= (D*dt)/(2*dx*dx);
k=1;



%defining IC & BC
IC= sin(k.*x);
u(length(t),length(x))=0;
u(1,:)=IC; %at time zero
g0=0;
gL=0;
%apply Ditchlet BC at all points in u...
u(:,1)=g0;
u(:,length(x))=gL;

% Given Closed-form Solution:
u_true(length(t),length(x))=0;
for time=1:1:length(t)
    for space=1:1:length(x)
        u_true(time,space)=exp(-D*k*k*t(time))*sin(k*x(space));
    end
end
%for Thomas Algorithm...
b= -lambda; %wont change-infact turned to zero on foward path
%For helping solve the tri-diagonal matrix
a(1,length(x))=0;
a(1,:)= 1+2*lambda; %This guy will be changed
c=-lambda; %will not change
y(length(x))=0;
A= 1+2*lambda;
C=lambda;

%Begin solving for the Diffusion equation
 for time=2:1:length(t)%for each time u 
 
        for space=2:1:(length(x)-1) %for each space u
            % first solve for u's known; will be y in the matrix
            y(space)= (u(time-1,space-1)*lambda)   +   ((-2*lambda+1)*u(time-1,space))  +  (u(time-1,space+1)*lambda);
        
        end
        for space=3:1:(length(x)-1) %Thomas foward path to alter a's and y's in matrix
            a(space)= a(space)-(b/a(space-1))*c;
            y(space)= y(space)-(b/a(space-1))*y(space-1);
        end
        
        u(time,length(x)-1)=y(length(x)-1)/a(length(x)-1); %first back-wards solve
        
        for space=(length(x)-2):-1:2 %solve rest of u's,completing Thomas back-wards 
         u(time,space)= (y(space)-c*u(time,space+1))/a(space);
        end
        %reset a values for tri-diagonal matrix for next time step, else solution will go to hell;
         a(1,length(x))=0;
         a(1,:)= 1+2*lambda;
        y(length(x))=0;
 end
 
    
       
        u_error=0;
        for time=2:1:length(t)-1
            for space=2:1:length(x)-1
                u_error= u_error +abs(u(time,space)-u_true(time,space));
            end
        end
        u_error= u_error/((length(t)-2)*(length(x)-2))
%plotting attemp
% uncomment  the for loop and the reset for an 'animated' plot
      for p=1:1:length(t)
%comment p=x for loop animation        
%         p=151;
          plot(x,u_true(p,:))
          hold on
          plot(x,u(p,:),'g--o')
          ylim([-1.2 1.2])
          title([' Problem one graph comparison at   ',num2str(t(p)), '  seconds'])
          xlabel('rod length')
          ylabel(' u magnitude')
          legend('closed-form','Numerical by Crank-Nicholson')
          pause(.01);
          clf('reset')
     
      end
        
%         mesh(u)
%         hold on
%         mesh(u_true)
%         title(' surface mesh of solutions')
%         ylabel('Time incremets')
%         xlabel(' Length of rod')
%         