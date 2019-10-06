%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Created by Aurelio Vasquez
%
%
% Computation of the first four moments of the NGARCH model
% using monte carlo simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
warning off MATLAB:nchoosek:LargeCoefficient;
 %option`s contract parameters
 dt=1/365;            %lenght of a time interval in year
 s0=50;               %initial stock price
 ks= [1.1; 1.0; 0.9;]; %strike to stock price ratio
 r=0.05*dt;           %risk free interest rate (per period)
 k=ks.*s0;            %transform ratio into strike prices
%[Binomial] = BinomialAmerican1(1, s0, s0, 0.05, 0, 270/365, .2, 100)
%GARCH process parameters

%[vol,fval]=fsolve(inline('ImpliedVol(1, s0, s0, 0.05, 0, 270/365, x0, 100,Binomial)','x0'),x0);
 b0=1e-5; b1=0.70; b2=0.1; lath=0.5; 

 T=30;
 scenarios=100000;
%initial variance
 lath2=(lath)^2; h1=(b0/(1-b2*(1+lath2)-b1)); %stationary variance

 epsilon=normrnd(0,1,T+1,scenarios);
 h=ones(T+1,scenarios)*h1;
 logRet=ones(T,scenarios);
 
 a=b1+b2.*(epsilon-lath).^2;
 b=h1*cumprod(a,1);

 inc=1:T+1;
 b0Mat=repmat(inc',1,scenarios);
 c=b0Mat*b0+b;
 for i=1:scenarios
 for j=2:31
    h(j,i)=b0+h(j-1,i)*(b1+b2*(epsilon(j-1,i)-lath)^2);   
    logRet(j-1,i)=r-.5*h(j,i)+h(j,i)^.5*epsilon(j,i);
 end
 end
 
 sret=sum(logRet,1)';
% hist(sret,scenarios^.5);

 
 mu1=mean(sret);
 mu2=std(sret);
 mu3=skewness(sret);
 mu4=kurtosis(sret);

 [f,r] = ksdensity(sret); 
 plot(r,f,'Color','r','LineWidth',2)
 hold on
 x=min(sret):.001:max(sret);
 y=normpdf(x,mu1,mu2);
 plot(x,y)
 hold off
 
 fprintf(' \n\n');             % line feed

  fprintf('First four moments : \n');
  fprintf('   %8.5f',mu1, mu2, mu3, mu4');

  fprintf(' \n\n');