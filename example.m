%
% An example on how to use the analytical approximation program
% for the GARCH option pricing model (European options). To use
% the program two files are required : gap.m and third_r2.dll 
%


 clear all;

%option`s contract parameters
 dt=1/365;            %lenght of a time interval in year
 s0=50;               %initial stock price
 ks=[0.9; 1.0; 1.1;]; %strike to stock price ratio
 r=0.05*dt;           %risk free interest rate (per period)
 k=ks.*s0;            %transform ratio into strike prices

%GARCH process parameters
 b0=1e-5; b1=0.7; b2=0.1; lath=1.0; 

%initial variance
 lath2=(lath)^2; h1=(b0/(1-b2*(1+lath2)-b1)); %stationary variance

%compute 30 and 90 days option prices 
 T=30;   [ca]=gap(b0,b1,b2,lath,h1,T,r,s0,k,dt,1);
 T=90;   [ca]=gap(b0,b1,b2,lath,h1,T,r,s0,k,dt,1);
 
