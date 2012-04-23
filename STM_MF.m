
%   STM_MF.m
%Forward Euler for MeanField e-e net
clear all; hold on; clc;
DT = 2;  %Time increment as fraction of time constant
Final_Time = 500;   %Final time value for calculation
Last = Final_Time/DT + 1;  %Last time step
Time = DT*[0:Last-1];  %Time vector
Tau = 20;  %Neural time constants in msec
X = zeros(1,Last);  %Vector to store response of Neuron #1
X(1)=30;
% Here's the Euler integration scheme
for T = 2:Last
 	X(T) = X(T-1)+ DT/Tau*(-X(T-1) + 100*(3*X(T-1))^2/(120^2 + (3*X(T-1))^2));  %Your Equation Here
end;
plot(Time,X);
xlabel('t'); ylabel('Firing rate');
axis([0 Final_Time 0 100]);

