
%   STM_MF.m
%Forward Euler for MeanField e-e net
clear all; hold on; clc;% clg;
global Tau = 20;  %Neural time constants in msec
global TauA = 2000;
global M = 100;
global sigma = 120;
global a = 3;
global I = 0;
global gamma=1;

function y=Snr(u)
  global M; global sigma;
  up=max(u,0);
  y=M*(up^2)/(sigma^2 + (up^2));
end

function y=tau_dr_dt(cur_R, cur_A)
  global a; global I;
  y= -cur_R + Snr(a*cur_R + I - cur_A);
end

function y=tau_da_dt(cur_R, cur_A)
  global gamma;
  y=-cur_A + gamma*cur_R;
end

function [Time, R_t, A_t] = getTimeSeries
  global Tau; global TauA;
  DT = 2;  %Time increment as fraction of time constant
  Final_Time = 15000;   %Final time value for calculation
  Last = Final_Time/DT + 1;  %Last time step
  Time = DT*[0:Last-1];  %Time vector

  R_t = zeros(1,Last);  %Vector to store response of Neuron #1
  R_t(1)=40; %ro
  A_t = zeros(1,Last);
  A_t(1) = 0;
  for T = 2:Last
    R_t(T) = R_t(T-1) + (DT/Tau)*tau_dr_dt(R_t(T-1), A_t(T-1));
    A_t(T) = A_t(T-1) + (DT/TauA)*tau_da_dt(R_t(T-1), A_t(T-1));
  end
end
TauA = 2000;
I = 30;
[Time, R_t, A_t] = getTimeSeries;
I = 50;
[Time, R_t2, A_t2] = getTimeSeries;
I = 70;
[Time, R_t3, A_t3] = getTimeSeries;

rVals = 0.0001:0.5:99.0001;
function x=solveQuad(b, c)
  a = 1;
  disc = b^2-4*a*c;
  x = (-b + sqrt(disc))/2;
end
function tIVals=steadyState(rVals)
  global a; global sigma; global M; global I;
  Iarr = zeros(1,size(rVals)(2)).+ I;
  nbVals = 2*a*(rVals);
  cVals = a^2*(rVals.^2) - (sigma^2)./((M./rVals).- 1);
  tIVals= arrayfun(@solveQuad, nbVals, cVals);
end
tIVals = steadyState(rVals);

function x=na(v,l)
  x = zeros(1,size(l)(2)).+v;
end

%plot(Time ,R_t,  ";r(t) with I=30;",Time, A_t,   ";A(t) with I=30;")
%xlabel('t'); ylabel('Firing rate');
%axis([0 Time(size(Time)(2)) 0 100]);

plot(Time ,R_t,  ";r(t) with I=30;",Time, A_t,   ";A(t) with I=30;",
     Time ,R_t2,";r(t) with I=50;", Time, A_t2, ";A(t) with I=50;",
     Time ,R_t3,";r(t) with I=70;", Time, A_t3, ";A(t) with I=70;");
xlabel('t'); ylabel('Firing rate');
axis([-60 Time(size(Time)(2)) 0 100]);

%plot(na(30,R_t) - A_t ,R_t,  ";r(t) with I=30;", na(30,R_t)- A_t, A_t,   ";A(t) with I=30;",
%     na(50,R_t)- A_t2 ,R_t2,";r(t) with I=50;", na(50,R_t)- A_t2, A_t2, ";A(t) with I=50;",
%     na(70,R_t)- A_t3 ,R_t3,";r(t) with I=70;", na(70,R_t)- A_t3, A_t3, ";A(t) with I=70;",
%     tIVals, rVals, ";Sytem steady state: dr/dt=0;");%,
%xlabel('I - A(t)'); ylabel('Firing rate');
%axis([-60 100 0 100]);

