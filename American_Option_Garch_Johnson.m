%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Created by Aurelio Vasquez
%
%
%  Function to price American Options using Johnson Curves
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Sn,Cn]=american_option_garch_Johnson(b0,b1,b2,la,h1,mat,r,s0,k,dt);
%
%  Computes the analytical approximation for the GARCH option pricing model
%  ------------------------------------------------------------------------
%
%  input   :  b0,b1,b2,la,h1  -- GARCH process parameters (b0 must be per period)
%             mat             -- maturity in number of periods
%             r               -- risk free rate per period
%             s0              -- initial stock price
%             k               -- vector of strike prices
%             
n=mat;   % Number of periods used.

%compute moment conditions
la2=la^2; la4=la^4; 
m1=b2*(1+la2)+b1; 
m2=(b2^2)*(3+6*la2+la4)+2*b2*b1*(1+la2)+(b1^2);
m3=(b2^3)*(15+45*la2+15*la4+(la^6))+3*b1*(b2^2)*(3+6*la2+la4)...
    +3*(b1^2)*b2*(1+la2)+(b1^3);

if m1>1 | m2>1 | m3>1
    error('GARCH parameter values are outside the admissible region');
end

if m1>0.98 | m2>0.98 | m3>0.98
    warning('GARCH parameter values are at the limit of the admissible region');
    warning('Results may be unreliable');
end


%compute the analytical moments of h 
[veh1]=cveh1(b0,b1,b2,la,h1,mat);
[veh2]=cveh2(b0,b1,b2,la,h1,mat);
[veh3]=cveh3(b0,b1,b2,la,h1,mat);

%compute the analytical approximate moments of cumulative returns 
[er1]=fm(b0,b1,b2,la,h1,mat,r,veh1);
[er2]=sm(b0,b1,b2,la,h1,mat,r,veh1,veh2);
[er3,vs1,vs2,vs3]=third_r2(b0,b1,b2,la,h1,mat,r,veh1,veh2,veh3);
[er4]=fom1(mat,r,vs1,vs2,vs3);

%standard deviation and variance of cumulative returns 
var=er2-(er1^2); stdv=sqrt(var); 

%analytical third moment of standardized return
s3=1/(stdv^3); mu3 = s3*er3 - 3*s3*er2*er1 + 3*s3*er1*(er1^2) - s3*(er1^3);

%analytical fourth moment of standardized returns
s4=1/(stdv^4); mu4 = s4*er4 - 4*s4*er3*er1 + 6*s4*er2*(er1^2) - 4*s4*(er1^3)*er1 + (s4*er1^4);

% Analytical approximation using the Edgeworth expansion. It uses the third and forth moments.
% Contruction of the binomial tree starts by considering an n-step binomial
% distribution denoted by Yj

r=r/dt;
T=mat/365;

j=0:1:n;
j1=j;
j1(1)=1;

%Construction of variables for our distrubution

b=(.5^n)*(factor(n)); %Factorial trick!!!
y=-3:6/n:3;
y=y*stdv+er1;

%Johnson Expansion

[gamma,delta,lambda,xi,itype,ifault] = johnson_pdf(er1,stdv,(mu3/stdv^3),(mu4/stdv^4));

x=xi+lambda.*exp((y-gamma)./delta);

if itype==1
    x=xi+lambda.*exp((y-gamma)./delta);
elseif itype==2
    x=xi+lambda.*sinh((y-gamma)./delta);
elseif itype==3
    x=(1+exp(-(y-gamma)/delta)).^-1;
end    

%Rescaling of probabilities to sum up to 1

P=b./sum(b); 
y=x;


%Mean and Variance
M=sum(P.*y);
V=sum(P.*(y-M).^2)^.5;

%Modified x with new probabilities
x=(y-M)/V;

% Computation of risk neutral probabilities
p=P./factor(n);

% Computation of stock price process
stdv=stdv/T^.5;

mu=r-1/T*log(sum(P.*exp(stdv.*x.*T^.5)));

S=s0.*exp(mu*T+stdv.*x.*T^.5);

%Backward recursion for stock price and option price. Since it's american,
%at each step I have to consider if the option is exercised or not.

C=max(k-S,0);
[Sn,Cn]=BinRec(p,S,C,k,r,T,n);



%print results

%   fprintf(' \n\n');
%   fprintf('Input values: \n');
%   fprintf('  r=%g,',r); 
%   fprintf('  stock price=%g,',s0); 
%   fprintf('  maturity=%g\n\n',mat);
% 
%   fprintf('GARCH param:  \n');
%   fprintf('  b0=%g,',b0);
%   fprintf('  b1=%g,',b1);
%   fprintf('  b2=%g,',b2);
%   fprintf('  lambda+theta=%g\n\n',la);
% 
%   fprintf('Strike prices: \n'); 
%   fprintf('   %8.5f',k);
% 
%   fprintf(' \n\n');             % line feed
% 
%   fprintf('Analytical approximation for European call prices : \n');
%   fprintf('   %8.5f',cap4');
% 
%   fprintf(' \n\n');

function [S,C]=BinRec(P,S,C,K,r,T,n);

% This function calculates the numerical derivatives of a function
% The function should return a Nx1 matrix and the parameter vector should
% be kx1. The returned value is a Nxk matrix

if size(P',1) ~= 1
    p=P(2:size(P',1))+P(1:size(P',1)-1);
    STemp = P.*S.*exp(-r*T/n);
    Snew=(STemp(2:size(STemp',1))+STemp(1:size(STemp',1)-1))./p;
    
    CStock=max(K-Snew,0);
    
    CTemp = P.*C.*exp(-r*T/n);
    Cnew=(CTemp(2:size(CTemp',1))+CTemp(1:size(CTemp',1)-1))./p;
    %Cnew=C(2:size(C',1))+C(1:size(C',1)-1);
    
    Cnew=max(Cnew,CStock);
    
    [S,C]=BinRec(p,Snew,Cnew,K,r,T,n);
    
    % else
    %     p=P(2:size(P',1))+P(1:size(P',1)-1);
    %     STemp = P.*S.*exp(-r*T/n);
    %     Snew=STemp(2:size(STemp',1))+Snew(1:size(STemp',1)-1)./P;
    %     CStock=max(S-K,0);
    %         
    %     CTemp = P.*C.*exp(-r*T/n);
    %     Cnew=(CTemp(2:size(CTemp',1))+CTemp(1:size(CTemp',1)-1))./P;
    %     Cnew=C(2:size(C',1))+C(1:size(C',1)-1);
    %     Cnew=max(Cnew,CStock);
    %     S=Snew;
    %     C=Cnew;
end


function [p]=factor(x);
% Computes a vector with the factorial equation
% (factorial(n)./(cumprod(j1).*fliplr(cumprod(j1))))
temp=0:1:x;
p(1)=1;
p(x+1)=1;
for i=2:ceil((x+1)/2)
    
    %   p(i)=prod(temp(1,x-i+3:x+1))/factorial(i-1);
    p(i)=nchoosek(x,i-1);
    p(x-i+2)=p(i);
    
end

function [eh]=cveh1(b0,b1,b2,lath,h1,mat);
%
%  Computes the conditional expectation of h(t) for t=1 to mat using the corrected
%  formula given in Computing American GARCH option prices using a Markov chain approximation
%
%  input   :  b0,b1,b2,lath,h1  -- GARCH process parameters (b0 must be per period)
%             mat               -- maturity in number of periods
%
%  output  :  eh                -- a vector of conditional expectation of h(t) for t=1 to T
%

%build an index vector to compute the powers  
lath2=(lath)^2;
index1=[0:(mat-1)]';
v=b1+b2*(1+lath2);
vt=v.^index1;

%compute the conditional expectation vector 
eh = h1.*vt + b0*((1-vt)./(1-v));



function [veh2]=cveh2(b0,b1,b2,lath,h1,mat);
%
%  Computes the second moment of h(t) for t=1 to t=T
%  -------------------------------------------------
%  
%  The variance is computed using the corrected formula given in 
%  American GARCH option pricing by a Markov Chain approximation 
%
%  input   :  b0,b1,b2,lath,h1 -- GARCH process parameters (b0 must be per period)
%             mat              -- maturity in number of periods
%


%values needed in the computations 
lath2=lath^2; b2_2=b2^2; b0_2=b0^2; 
vv=b2*(1+lath2)+b1; 
uu=(vv^2) + 2*(1+2*lath2)*b2_2;

%build an index and compute powers of vv and uu 
index1=[0:(mat-1)]';

ut=uu.^index1;
vt=vv.^index1;

ra1=(vv.*(ut-vt)) ./ (uu-vv);
ra2=(1-ut) ./ (1-uu);
ra3=vv / (uu-vv);
ra4=(1-ut) ./ (1-uu);
ra5=(1-vt) ./ (1-vv);

part1 = (h1^2).*ut;
part2 = 2*b0*h1 .* ra1;
part3 = b0_2 .* (ra2 + 2*ra3 .* (ra4-ra5));

veh2 = part1 + part2 + part3;



function [veh3]=cveh3(b0,b1,b2,lath,h1,mat);
%
%  Computes the third moment of h(t) for t=1 to t=T
%  -------------------------------------------------
%  
%  The variance is computed using the corrected formula given in 
%  American GARCH option pricing by a Markov Chain approximation 
%
%  input   :  b0,b1,b2,lath,h1 -- GARCH process parameters (b0 must be per period)
%             mat              -- maturity in number of periods
%


%compute values for mu1, mu2 and mu3
lath2=lath^2; lath4=lath^4; lath6=lath^6;
m1=b2*(1+lath2)+b1; 
m2=(b2^2)*(3+6*lath2+lath4)+2*b2*b1*(1+lath2)+(b1^2);
m3=(b2^3)*(15+45*lath2+15*lath4+lath6)+3*b1*(b2^2)*(3+6*lath2+lath4)...
    +3*(b1^2)*b2*(1+lath2)+(b1^3);

%build an index and compute powers for the mu's
index=[0:(mat-1)]';
m1t=m1.^index; m2t=m2.^index; m3t=m3.^index;

%compute the ratios needed in the computations
r1=(1-m3t)./(1-m3);
r2=m2./(1-m2);
r3=(m2t-m3t)./(m2-m3);
r4=m1./(1-m1);
r5=(m1t-m3t)./(m1-m3);
r6=m2./(m1-m2);

%build the different parts of the equation
p1=r1;
p2=3.*(r2.*r1-r2.*r3);
p3=3.*(r4.*r1-r4.*r5);
p4=6.*(r4.*r2.*r1-r4.*r2.*r3-r4.*r6.*r5+r4.*r6.*r3);
part1=(b0^3).*(p1+p2+p3+p4);

p1=r5;
p2=2.*(r6.*r5-r6.*r3);
part2=3.*(b0^2).*m1.*(p1+p2).*h1;

part3=3.*b0.*m2.*r3.*(h1^2) + m3t.*(h1^3);

veh3 = part1 + part2 + part3;



function [er1]=fm(b0,b1,b2,la,h1,mat,r,veh1);
%
%  Computes first moment of return
%  -------------------------------
%
%  input   :  b0,b1,b2,la,h1 -- GARCH process parameters
%             mat            -- maturity in number of periods
%             r              -- risk free rate per period
%             veh1           -- expected values of h(t) for 
%                               t=1 to mat
%

%compute expected value
sveh1 = sum(veh1);  
er1 = -0.5*sveh1 + mat*r;



function [er2]=sm(b0,b1,b2,lath,h1,t,r,veh1,veh2);
%
%  Computes the second moment of return
%  ------------------------------------
%
%  input   :  b0,b1,b2,lath,h1  -- GARCH process parameters (b0 must be per period)
%             t                 -- maturity in number of periods
%             r                 -- risk free rate (per period)
%             veh1              -- vector of expected values for h(t)^1 for t=1 to mat 
%             veh2              -- vector of expected values for h(t)^2 for t=1 to mat
%


%compute summations

%values required in the computations 
lath2=lath^2; y=b2*(1+lath2)+b1; 
index=[1:t]'; tmi=t-index; yt=y.^tmi; omy=1-y; omyt=1-yt;

%compute first part
part1=(t^2)*(r^2);

%second part
part2=t*r*sum(veh1);

%third part
r1=tmi./(omy); r2=y./omy; r3=omyt./omy; 
p1=b0*sum( (r1-r2.*r3).*veh1 );
p2=sum( y.*r3.*veh2 );
p3=sum(veh2);

sd1=2*(p1+p2)+p3;

part3=0.25*sd1;

%fourth part
part4=sum(veh1);

%fifth part
[veh32]=nim(b0,b1,b2,lath,h1,t,veh1,veh2,3);
part5=-2*lath*b2*sum(  r3 .* veh32 );

%compute second moment 
er2=part1 - part2 + part3 + part4 - part5;



function [ehni]=nim(b0,b1,b2,lath,h1,mat,eh1,eh2,n);
%
%  Computes non integer moments of h(t) using a an approximation based
%  on a Taylor series truncated after the second term
%
%  input   :  b0,b1,b2,lath,h1 -- GARCH process parameters (b0 must be per period)
%             mat              -- maturity in number of periods
%             eh1              -- a vector of expected values for h^1
%             eh2              -- a vector of expected values for h^2
%             n                -- integer value in n/2 (the non-integer power)
%
%  output  :  ehni             -- non-integer moments of h(t) for t=1 to T

%compute approximation for the non-integer moment
po=n/2; pou=po-1; pod=pou-1; 
ehni = (eh1 .^po) + (.5)*po*pou*(eh1 .^pod) .*(eh2-(eh1 .^2));






function [er4]=fom1(t,r,vs1,vs2,vs3);
%
%  Computes the fourth conditional moment of return (with some terms left out)
%  ---------------------------------------------------------------------------
%
%  input   :  t              -- maturity in number of periods
%             r              -- risk free interest rate (per period)
%             vs1,vs2,vs3    -- vector of values extracted from the c-proc third.c
%
%  output  :  er4            -- Fourth moment of stock return
%
%  globals :  none
%


%extract values that are available from the third moment computation
sum1=vs1(1,1); sum2=vs1(2,1); sum3=vs1(3,1); sum4=vs1(4,1); sum5=vs1(5,1); 
sum6=vs1(6,1); sum7=vs1(7,1); sum8=vs1(8,1); sum9=vs1(9,1); sm10=vs1(10,1); 
sm11=vs1(11,1); sm12=vs1(12,1); q5_8=vs1(13,1);

sa=vs2(1,1); sb=vs2(2,1); sc=vs2(3,1); sd=vs2(4,1); 
se=vs2(5,1); sf=vs2(6,1); sg=vs2(7,1); 

sumh=vs3(1,1); sumh2=vs3(2,1); sumh32=vs3(3,1); sumh3=vs3(4,1);

%compute quantities required in the fourth moment computation
soee=12*sum8+6*sum6+3*sumh2;
sq4 = 2*sm11 +  2*sum9;
sq5 = 3*sum7 + 3*sm12 + 3*sum3 + 3*sm10 + q5_8;

%compute the different parts of the fourth moment 
p1 = (t^4)*(r^4);
p5 = -2*(t^3)*(r^3)*sumh;
p6 = ((3/2))*(t^2)*(r^2)*sc;
p9 = 6*(t^2)*(r^2)*se;
p10 = 4*t*r*sb;
p12 = -6*(t^2)*(r^2)*sg;
p13 = 3*t*r*sd;
p14 = -6*t*r*sf;
p15 = soee;
p16 = -t*r*sa;

%compute the fourth moment
er4 = p1+p5+p6+p9+p10+p12+p14+p15+p16 + (3/2)*sq4 - 2*sq5;





function [pro]=cnor(x);
%
%  Compute cumulative normal distribution function
%  -----------------------------------------------
%
%  ref: Hull, p. 226
%

%initialize the prob. matrix to zero
pro=zeros(size(x));  

%initialise values needed in computations 
a1= 0.319381530; a2=-0.356563782; a3= 1.781477937; a4=-1.821255978; a5=1.330274429; ga=0.2316419;

%cases for which x >= 0
ind=find(x>=0);
k  = 1 ./ (1+(ga.*x(ind)));
np = (1/(sqrt(2*pi))) .* exp(-((x(ind).^2)/2));  
pro(ind)=1 - np.*(a1.*k+a2.*(k.^2)+a3.*(k.^3)+a4.*(k.^4)+a5.*(k.^5));

%cases for which x < 0
ind=find(x<0); 
mx=-x(ind);  
k = 1 ./ (1+(ga.*mx));
np = (1/(sqrt(2*pi))) .* exp(-((mx.^2)/2)); 
umpro=1 - np.*(a1.*k+a2.*(k.^2)+a3.*(k.^3)+a4.*(k.^4)+a5.*(k.^5));
promx=umpro; 
pro(ind)=1-promx;



function [gamma,delta,lambda,xi] = johnson(mu1,mu2,mu3,mu4);

b1=(mu3/mu2^3)^2;
b2=mu4/mu2^4;
b3=b2-3;
w = roots([1 1 2 -4-b1]);

b22=w.^4+2*w.^3+3*w.^2-3;
for i=1:3
    if imag(w(i,1))==0
        if b22(i,1)<b2 
            method=1;
            w1=w(i,1);
            w=w1;
            
        else 
            w1=w(i,1);
            w=w1;
            method=2;
            if imag(((2*b2-2.8*b1-2)^.5-1)^.5)~=0
                method=1;
            end
        end 
    end
end

if method==1 
    delta=(log(w))^-.5;
    gamma=.5*(delta*log(w*(w-1)/mu2));
    xi=sign(mu3)*mu1-exp((1/(2*delta)-gamma));
    lambda=sign(mu3);
elseif method==2
    if b1==0
        delta=(log(w))^-.5;
        gamma=0;
    else
        w=((2*b2-2.8*b1-2)^.5-1)^.5;
        z=100;
        while abs(b1 - z) > 0.01
            if(abs((b1)^.5)>0.01)
                w1=w+1;
                wm1 = w - 1;
                z = w1 * b3;
                v = w * (6 + w * (3 + w));
                a = 8 * (wm1 * (3 + w * (7 + v)) - z);
                b = 16 * (wm1 * (6 + v) - b3);
                y = ((a * a - 2 * b * (wm1 * (3 + w * (9 + w * (10 + v))) - 2 * w1 * z)) ^.5- a) / b;
                z = y * wm1 * (4 * (w + 2) * y + 3 * w1 * w1) ^ 2 /  (2 * (2 * y + w1) ^ 3);
                v = w * w;
                w = (1 - 2 * (1.5 - b2 + (b1 *  (b2 - 1.5 - v * (1 + 0.5 * v))) / z))^.5;
                w = (w-1);
            end
        end
        y = y / w;
        y = log(sqrt(y) + sqrt(y + 1));
        if (rb1 > 0) 
            y = -y;
        end
        x = sqrt(1 / log(w));
        delta = x;
        gamma = y * x;
        y = exp(y);
        z = y * y;
        x = mu2 / sqrt(0.5 * (w - 1) * (0.5 * w *  (z + 1 / z) + 1));
        lambda = x;
        xi = (0.5 * sqrt(w) * (y - 1 / y)) * x + mu1;
        
        
    end
end 
