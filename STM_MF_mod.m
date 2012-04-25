
%   STM_MF.m
%Forward Euler for MeanField e-e net
clear all; hold on; clc;% clg;
global Tau = 20;  %Neural time constants in msec
global TauA = 2000;
global M = 100;
global sigma = 120;
global a = 3;
global I = 0;
global gamma=2;

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

function x=na(v,l)
  x = zeros(1,size(l)(2)).+v;
end

function first()
  global gamma;
  gamma = 2;
  [Time, R_t, A_t] = getTimeSeries;
  gamma = 1;
  [Time, R_t2, A_t2] = getTimeSeries;

  plot(Time ,R_t, ";r(t) with gamma=2;", Time, A_t,  ";A(t) with gamma=2;",
       Time ,R_t2,";r(t) with gamma=1;", Time, A_t2, ";A(t) with gamma=1;")
  xlabel('t'); ylabel('Firing rate');
  axis([0 2000 0 150]);
  print("1.png")
  gamma = 2; %Set control back
end
function second()
  global TauA;
  TauA = 2000;
  [Time, R_t, A_t] = getTimeSeries;
  TauA = 8000;
  [Time, R_t2, A_t2] = getTimeSeries;
  TauA = 16000;
  [Time, R_t3, A_t3] = getTimeSeries;
  TauA = 2000;

  plot(Time ,R_t, ";r(t) with TauA=2000;", Time, A_t,  ";A(t) with TauA=2000;",
       Time ,R_t2,";r(t) with TauA=8000;", Time, A_t2, ";A(t) with TauA=8000;",
       Time ,R_t3,";r(t) with TauA=16000;",Time, A_t3, ";A(t) with TauA=16000;")
  xlabel('t'); ylabel('Firing rate');
  axis([0 Time(size(Time)(2)) 0 100]);
  print("2.png")
end

function third()
    global TauA;
  TauA = 2000;
  [Time, R_t, A_t] = getTimeSeries;
  TauA = 8000;
  [Time, R_t2, A_t2] = getTimeSeries;
  TauA = 16000;
  [Time, R_t3, A_t3] = getTimeSeries;
  TauA = 2000;

  rVals = 0.0001:0.5:99.0001;
  tIVals = steadyState(rVals);

  plot(-A_t ,R_t, ";r(t) with TauA=2000;",
       -A_t ,R_t2,";r(t) with TauA=8000;", 
       -A_t ,R_t3,";r(t) with TauA=16000;", -A_t, A_t,  ";A(t);",
       tIVals, rVals, ";System steady state: d/dt=0;")
  xlabel('I - A(t)'); ylabel('Firing rate');
  axis([-60 20 0 80]);
  print("3.png")
end

function fourth()
  global I;

  rVals = 0.0001:0.5:99.0001;

  I = 30;
  [Time, R_t, A_t] = getTimeSeries;
  tIVals1 = steadyState(rVals);
  I = 50;
  [Time, R_t2, A_t2] = getTimeSeries;
  tIVals2 = steadyState(rVals);
  I = 70;
  [Time, R_t3, A_t3] = getTimeSeries;
  tIVals3 = steadyState(rVals);
  I = 0;

  plot(Time ,R_t,  ";r(t) with I=30;",Time, A_t,   ";A(t) with I=30;")
  xlabel('t'); ylabel('Firing rate');
  axis([0 Time(size(Time)(2)) 0 120]);
  print("4ta.png")
  clf

  plot(Time ,R_t2,";r(t) with I=50;", Time, A_t2, ";A(t) with I=50;")
  xlabel('t'); ylabel('Firing rate');
  axis([0 Time(size(Time)(2)) 0 120]);
  print("4tb.png")
  clf
  plot(Time ,R_t3,";r(t) with I=70;", Time, A_t3, ";A(t) with I=70;")
  xlabel('t'); ylabel('Firing rate');
  axis([0 Time(size(Time)(2)) 0 120]);
  print("4tc.png")
  clf

  plot(na(30,R_t) - A_t ,R_t,  ";r(t) with I=30;", na(30,R_t)- A_t, A_t,   ";A(t) with I=30;",
       tIVals1, rVals, ";System steady state for I = 30: d/dt=0;")
  xlabel('I - A(t)'); ylabel('Firing rate');
  axis([-60 100 0 100]);
  print("4pa.png")
  clf

  plot(na(50,R_t)- A_t2 ,R_t2,";r(t) with I=50;", na(50,R_t)- A_t2, A_t2, ";A(t) with I=50;",
       tIVals2, rVals, ";System steady state for I = 50: d/dt=0;")
  xlabel('I - A(t)'); ylabel('Firing rate');
  axis([-60 100 0 100]);
  print("4pb.png")
  clf
  plot(na(70,R_t)- A_t3 ,R_t3,";r(t) with I=70;", na(70,R_t)- A_t3, A_t3, ";A(t) with I=70;",
       tIVals3, rVals, ";System steady state for I = 70: d/dt=0;")
  xlabel('I - A(t)'); ylabel('Firing rate');
  axis([-60 100 0 100]);
  print("4pc.png")
  clf
end

function fifth(low,high)
  global I;
  I = low;
  [Time, R_t1, A_t1] = getTimeSeries;
  I = high;
  [Time, R_t2, A_t2] = getTimeSeries;
  I = 0;

  clf
  plot(Time, R_t1,strcat(";r(t) with I=",num2str(low),";"),
       Time, A_t1,strcat(";A(t) with I=",num2str(low),";"),
       Time ,R_t2,strcat(";r(t) with I=",num2str(high),";"),
       Time, A_t2,strcat(";A(t) with I=",num2str(high),";"))
  xlabel('t'); ylabel('Firing rate');
  axis([0 Time(size(Time)(2)) 0 120]);
  %print("5.png")
end
  
%TauA = 2000;

%
%plot(Time ,R_t,  ";r(t) with I=30;",Time, A_t,   ";A(t) with I=30;")
%xlabel('t'); ylabel('Firing rate');
%axis([0 Time(size(Time)(2)) 0 100]);


%plot(na(30,R_t) - A_t ,R_t,  ";r(t) with I=30;", na(30,R_t)- A_t, A_t,   ";A(t) with I=30;",
%     na(50,R_t)- A_t2 ,R_t2,";r(t) with I=50;", na(50,R_t)- A_t2, A_t2, ";A(t) with I=50;",
%     na(70,R_t)- A_t3 ,R_t3,";r(t) with I=70;", na(70,R_t)- A_t3, A_t3, ";A(t) with I=70;",
%     tIVals, rVals, ";Sytem steady state: dr/dt=0;");%,
%xlabel('I - A(t)'); ylabel('Firing rate');
%axis([-60 100 0 100]);

