clear all
clc
% Solution for HWk 8 H.E/ Diffusion equation
% part two only
%This is a MATLAB FILE
L=pi;
D=.1;
T=10;
F=0;
time =0;
space=0;

%will divide interval into steps equal  INNER steps
steps=40;
dx=L/(steps+1);
x=0:dx:pi;
%space interval into Twenty steps
time_steps=8000;
dt=T/time_steps;
t=0:dt:T;

%defining IC & BC
k=1;
IC= 0;
u(length(t),length(x))=0;
u(1,:)=IC; %at time zero
%For BC
% Here w is omega
DT=(1/dt);
%one omega or the other; NOT BOTH!!!!!!
%w=.1/dt;
w=3.5;
% g0=sin(w.*t);
% gL=sin(w.*t)*cos(k*L);
%apply Ditchlet BC at all points in u...
for l=2:1:length(t)
    
    u(l,1)=sin(w*t(l));     
    u(l,length(x))=sin(w*t(l))*cos(k*L)  ;
        

end



% Given Closed-form Solution:
u_true(length(t),length(x))=0;
for time=1:1:length(t)
    for space=1:1:length(x)
        u_true(time,space)=sin(w*t(time))*(cos(k*x(space)));
    end
end
%for Thomas Algorithm...
% defining 'lambda' to help solve crank Nicholson along w/ Thomas
lambda= (D*dt)/(2*dx*dx);
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
            %solving for given F(x) in the diffusion equation
              if space==2 % left end of equation...% first solve for u's known; will be y in the matrix
                  y(space)= (u(time-1,space-1)*lambda) +lambda*u(time,space-1)  +   ((-2*lambda+1)*u(time-1,space))  +  (u(time-1,space+1)*lambda) + dt*(   w*cos(w*t(time))+((D*k*k))*sin(w*t(time))  )*( cos(k*x(space)));   
              elseif space==length(x)-1  %last inner node must include the Bc to calculate this one
                   y(space)= (u(time-1,space-1)*lambda) +lambda*u(time,space+1)  +   ((-2*lambda+1)*u(time-1,space))  +  (u(time-1,space+1)*lambda) + dt*(   w*cos(w*t(time))+((D*k*k))*sin(w*t(time))  )*( cos(k*x(space)));
              else
                   y(space)= (u(time-1,space-1)*lambda)   +   ((-2*lambda+1)*u(time-1,space))  +  (u(time-1,space+1)*lambda) + dt*(   w*cos(w*t(time))+((D*k*k))*sin(w*t(time))  )*( cos(k*x(space)));   
              end
        end
        for space=3:1:(length(x)-1) %Thomas foward path to alter b's and y's in matrix
            a(space)= a(space)-(b/a(space-1))*c;
            y(space)= y(space)-(b/a(space-1))*y(space-1);
        end
        
        u(time,length(x)-1)=y(length(x)-1)/a(length(x)-1); %first back-wards solve i.e. the last one...
        
        for space=(length(x)-2):-1:2 %solve rest of u's,completing Thomas back-wards 
         u(time,space)= (y(space)-c*u(time,space+1))/a(space);
        end
        
        %reset a values for tri-diagonal matrix for next time step, else solution will go to hell;
        a(1,length(x))=0;
        a(1,:)= 1+2*lambda;
        y(length(x))=0;
 end
 
        u_error_matrix= abs(u-u_true);
        u_error=0;
        for time=2:1:length(t)-1
            for space=2:1:length(x)-1
                u_error= u_error +abs(u(time,space)-u_true(time,space));
            end
        end
        u_error= u_error/((length(t)-2)*(length(x)-2))
      %plotting attemp
      
%       for p=1:1:length(t)
%           p=2001;  %will be for time; make sure is less than length of t-vector    
%           plot(x,u_true(p,:))
%           hold on
%           plot(x,u(p,:),'g--o')
%           ylim([-1.2 1.2])
%           legend('closed-form','Numerical by Crank-Nicholson');
%           title([' Problem one graph comparison at   ',num2str(t(p)), '  seconds'])
%           xlabel('rod length')
%           ylabel(' u magnitude')
%          
%         pause(.0000001);
%         clf('reset')
%       end
        mesh(u)
        hold on
        mesh(u_true)
        title(' surface mesh of solutions')
        ylabel('Time incremets')
        xlabel(' Length of rod')
        
        