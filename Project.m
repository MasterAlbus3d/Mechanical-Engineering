% MATLAB Code for problem Crank Nickelson method
% Finite volume method for 2D transient conduction heat transfer in a bus roof
clear; clc; close all;

%geometry definition

L = 0.04; %thickness
NL = 15;
tmax = 60*60*24*4;
Nt = 60*15;
dL = L /(NL-1);
dt = tmax / (Nt);
%external conditions 
h1 = 37;
h2 = 8;
Ta2 = 30; 
%materials

% % 1- polysterene
% k = 0.125;
% rho = 1050;
% c = 1720;
% alpha = k/(rho*c);

% 2- Rigid Polymer foam
% rho =165;
% k =0.03	;
% c =1190;
% alpha = k/(rho*c);

% 3- Ceramic Foam
% rho =670;
% k =0.7;
% c =960;
% alpha = k/(rho*c);

% 4- Wool
rho =1320;
k =0.3;
c =1370;
alpha = k/(rho*c);

% 5- polyethylene
% rho =960;
% k =0.42;
% c =1850;
% alpha = k/(rho*c);

%Initial and boundary coniditons

T = zeros(NL,Nt);
TOL = 1e-7;
error = 1; iterations = 0;

T(1,1:Nt) = 30;  %Bottom B.C.
T(NL,1:Nt) = 45;  %Top B.C.

while iterations < 100 
        iterations = iterations+1;
          Told = T;
          
          for i = 2:NL-1    %interal nodes
              for j= 1:Nt-1
                  Ta1 = 32+12*sin((pi/12)*((j/4)-9));
                  Qrad = 950*cos((pi/12)*(j/4)-12);
    
                  T(1,j+1) = (2*h2*alpha*dt*dL*Ta2 - (2*k*alpha*dt + h2*alpha*dt*dL - 2*k*dL^2)*T(1,j)+ 2*k*alpha*dt*(T(2,j)+T(2,j+1)))/(h2*alpha*dt*dL + 2*k*(alpha*dt+ dL^2));
                    T(i,j+1) = (1/(2*alpha*dt+dL^2))*(alpha*dt*T(i-1,j)+ alpha*dt*T(i-1,j+1) - 2*alpha*dt*T(i,j) + alpha*dt*T(i+1,j) + alpha*dt*T(i+1,j+1) + dL^2*T(i,j));
                    T(NL,j+1) = (2*alpha*dt*(h1*dL*Ta1+k*(Qrad*dL+T(NL-1,j)+T(NL-1,j+1)))-(2*k*alpha*dt+h1*alpha*dt*dL-2*k*dL^2)*T(NL,j))/(h1*alpha*dt*dL+2*k*(alpha*dt+dL^2));
              end
          end
          
          error = max(max(abs(Told-T)));
end

subplot(2,1,1),contour(T),
title('Temperature (Transient)'),xlabel('time (s)'),ylabel('thickness'),colorbar
subplot(2,1,2),pcolor(T),shading interp,
title('Temperature (Transient)'),xlabel('time (s)'),ylabel('thickness'),colorbar
