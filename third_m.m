function [er3m,vs1m,vs2m,vs3m]=third_m(b0,b1,b2,la,h1,T,r,veh1,veh2,veh3);
%
% Computes the third moment and return quantities required to compute
% an approximate fourth moment
%
%  inputs  :  b0,b1,b2,la,h1    -- GARCH process parameters (b0 must be per period)
%             T                 -- maturity in number of periods
%             r                 -- risk free rate (per period)
%             veh1              -- vector of expected values for h(t)^1 for t=1 to T 
%             veh2              -- vector of expected values for h(t)^2 for t=1 to T
%             veh3              -- vector of expected values for h(t)^3 for t=1 to T
%
%  outputs :  er3m              -- third moment of return
%
%             vs1m              -- 13x1 vector. The elements of this vector are 
%                                  sum1 to sum12 and SQ5_8 (see description bellow)
%
%             vs2m              -- 7x1 vector of quantities required in the computation
%                                  of the fourth moment. These quantities are computed
%                                  using the elements of vs1m. The elements are (using
%                                  the notation of the published paper) :
%                                  ST1, ST2, SD1, ST3, SD2, ST4, SD3
%
%             vs3m              -- 4x1 vector. 
%                                  vs3m(1) is sum over t of E[h(t)]
%                                  vs3m(2) is sum over t of E[h(t)^2]
%                                  vs3m(3) is sum over t of E[h(t)^(3/2)]
%                                  vs4m(4) is sum over t of E[h(t)^3]
%
%

%quantities required in the computation of the single and double sums
 la2=la^2; la4=la^4; 

 mu1=b2*(1+la2)+b1;                                                                           %equation (41)
 mu2=(b2^2)*(3+6*la2+la4)+2*b2*b1*(1+la2)+(b1^2);                                             %equation (42)
 mu3=(b2^3)*(15+45*la2+15*la4+(la^6))+3*b1*(b2^2)*(3+6*la2+la4)+3*(b1^2)*b2*(1+la2)+(b1^3);   %equation (43)

 v1=-2*b2*la;                                                              %equation (44)
 v2=-4*b2*la*(b1+3*b2+b2*la2);                                             %equation (45)
 v3=-6*b2*b2*b2*la*(15+10*la2+la4) - 12*b1*b2*b2*la*(3+la2)-6*b1*b1*b2*la; %equation (46)
 k1=b1+b2*la2+3*b2;                                                        %equation (47)
 k2=b1*b1+2*b1*b2*(3+la2)+b2*b2*(15+18*la2+la4);                           %equation (48)
 e1=-6*b2*la;                                                              %equation (49)
 
%  fprintf(' mu1                   : %12.6f \n',mu1);
%  fprintf(' mu2                   : %12.6f \n',mu2);
%  fprintf(' mu3                   : %12.6f \n',mu3);
%  fprintf(' nu1                   : %12.6f \n',v1);
%  fprintf(' nu2                   : %12.6f \n',v2);
%  fprintf(' nu3                   : %12.6f \n',v3);
%  fprintf(' zeta1                 : %12.6f \n',k1);
%  fprintf(' zeta2                 : %12.6f \n',k2);
%  fprintf(' xi1                   : %12.6f \n',e1);
 
%compute the simple sums
%=======================

 %Note : Many of the double sums have been simplified into double sums. The algebra of the simplifications
 %       is explained in the document : Further Details on the implementation of the formulas in An Analytical 
 %       Approximation for the GARCH Option Pricing model
 %


%index used in the computation of the simple sums 
 i=(1:T)';  Tmip1=T-i+1;  Tmi=T-i; 

%constants used in the computation of the sums
 %c1=-(-T+i+mu1.*T-i.*mu1+mu1-mu1.^Tmip1)./((-1+mu1).^2);
 c1=(Tmi/(1-mu1))-(mu1/(1-mu1))*((1-mu1.^Tmi)/(1-mu1));
 c2=mu1.*((1-mu1.^Tmi)./(1-mu1));
 c3=(1-mu1.^Tmi)./(1-mu1);
 c4=(1./(1-mu2)).*(1./(1-mu1))+((mu2.^Tmip1)./((mu1-mu2)*(1-mu2)))-((mu1.^Tmip1)./((mu1-mu2)*(1-mu1)));
 c5=(1-mu2.^Tmi)./(1-mu2);
 c6=((1+mu1)/(1-mu1)) .* ((Tmi)./(1-mu2)) - 2.*((mu1^2)./(1-mu1)).*(1./(mu1-mu2)).*( (1-mu1.^Tmi)./(1-mu1) );
 c7=( ((2*mu1)/(1-mu1)) .* (mu2./(mu1-mu2)) - ((1+mu1)./(1-mu1)).*(mu2./(1-mu2)) ) .* ( (1-mu2.^Tmi)./(1-mu2) );

%sum1 : Double sum of h(i)h(i+j)
 sum1_1=sum(veh1.*c1); sum1_2=sum(veh2.*c2); sum1=b0.*sum1_1 + sum1_2;

%sum2 : Double sum of h(i)e(i)h(i+j)
 sum2=v1.*sum(veh2.*c3);
 
%sum3 : Double sum of h(i)^(3/2)e(i)h(i+j)
 [veh52]=nim(b0,b1,b2,la,h1,T,veh1,veh2,5);  sum3=v1.*sum(veh52.*c3);

%sum4 : Double sum of h(i)^(1/2)e(i)h(i+j)
 [veh32]=nim(b0,b1,b2,la,h1,T,veh1,veh2,3);  sum4=v1.*sum(veh32.*c3);

%sum5 : Double sum of h(i)^(1/2)h(i+j)^2
 sum5_1=sum(veh2.*c1); sum5_2=sum(veh3.*c2);  sum5=b0.*sum5_1+sum5_2;

%sum6 : Double sum of h(i)e(i)^2 h(i+j)
 sum6_1=sum(veh1.*c1); sum6_2=sum(veh2.*c3);  sum6=b0.*sum6_1+k1.*sum6_2;

%sum9 : Double sum of h(i)h(i+j)^2
 sum9_1=sum(veh1.*(c6+c7)); sum9_2=sum(veh2.*c4); sum9_3=sum(veh3.*mu2.*c5); sum9=b0^2*sum9_1+2*b0*mu1*sum9_2+sum9_3;

%sum10 : Double sum of h(i)^(1/2)e(i)h(i+j)^2
 sum10_1=sum(veh32.*c4); sum10_2=sum(veh52.*c5); sum10=2*b0.*v1.*sum10_1+v2.*sum10_2;

%last term of S(Q5)
 sq5_8=sum(veh52.*e1.*c3);


%compute the double sums
%-----------------------

 %Note : Many of the triple sums have been simplified into double sums. The algebra of the simplifications
 %       is explained in the document : Further Details on the implementation of the formulas in An Analytical 
 %       Approximation for the GARCH Option Pricing model.
 %


%initialize quantities required in the computation of the double sums

  %sum7 is the triple sum of h(i)h(i+j)^1/2e(i+j)h(i+j+k)
   sum7_1=0; sum7_2=0; sum7_3=0;                %sum7_1 is the first part of sum7, sum7_2 is the second, etc.

  %sum8 is the triple sum of h(i)^1/2e(i)h(i+j)^1/2e(i+j)h(i+j+k)
   sum8_1=0; sum8_2=0; sum8_3=0;                %sum8_1 is the first part of sum8, sum8_2 is the second, etc.

  %sum11 is the triple sum of h(i)h(i+j)h(i+j+k)
   sum11_1=0; sum11_2=0; sum11_3=0; sum11_4=0; sum11_5=0;  %sum11_1 is the first part of sum11, sum11_2 is the second, etc.

  %sum12 is the triple sum of h(i)^1/2e(i)h(i+j)h(i+j+k)
   sum12_1=0; sum12_2=0; sum12_3=0;
 

%compute the double sums
 for i=1:T

     %index for the double sums
      j=(1:T-i)'; Tmimj=T-i-j; ipj=i+j;

     %compute constant required in the sums
      c1=(mu1.^(Tmimj)-1)./(mu1-1);
      c2=(1-mu1.^j)./(1-mu1);
      c3=((1+mu1)/(1-mu1)) .* ((1-mu2.^j)./(1-mu2)) - 2.*(mu1./(1-mu1)).*((mu1.^j-mu2.^j)./(mu1-mu2));
      c4=(mu1.^j-mu2.^j)./(mu1-mu2);
      c5=(Tmimj./(1-mu1))-(mu1./(1-mu1)).*((1-mu1.^Tmimj)./(1-mu1));
      c6=(1-mu1.^Tmimj)./(1-mu1);
      

     %sum7
      sum7_1=sum7_1+sum( (1/8) .* c6 .* veh1(ipj).^(3/2) .* veh1(i) );

      sum7_2=sum7_2+sum( (3/4) .* b0 .* c2 .* c6 .* veh1(ipj).^(1/2) .* veh1(i) + ...
                         (3/4) .* (mu1.^j) .* c6 .* veh1(ipj).^(1/2) .* veh2(i)   );

      sum7_3=sum7_3+sum( (3/8) .* c3 .* c6 .* veh1(ipj).^(-1/2) .* veh1(i) .* (b0^2) + ...
                         (3/8) .* mu1 .* c4 .* c6 .* veh1(ipj).^(-1/2) .* 2 .* b0 .* veh2(i) + ...
                         (3/8) .* mu2.^j .* c6 .* veh1(ipj).^(-1/2) .* veh3(i) );

     %sum8
      sum8_1=sum8_1+sum( (3/4) .* c6 .* c4 .* veh1(ipj).^(-1/2) .* veh32(i) );
      sum8_2=sum8_2+sum( (3/8) .* mu2.^(j-1) .* c6 .* veh1(ipj).^(-1/2) .* veh52(i) );
      sum8_3=sum8_3+sum( (3/4) .* mu1.^(j-1) .* c6 .* veh1(ipj).^(1/2) .* veh32(i) );


     %sum11
      sum11_1=sum11_1+sum( veh1(i).*b0.*c2.*c5 );
      sum11_2=sum11_2+sum( veh2(i).*mu1.^j.*c5 );
      sum11_3=sum11_3+sum( veh1(i).*b0^2.*c3.*mu1.*c6 );
      sum11_4=sum11_4+sum( veh2(i).*2.*b0.*mu1.*c4.*mu1.*c6);
      sum11_5=sum11_5+sum( veh3(i).*mu2.^j.*mu1.*c6);

     %sum12
      sum12_1=sum12_1+sum( veh32(i).*mu1.^(j-1).*c5 );
      sum12_2=sum12_2+sum( veh32(i).*c4.*mu1.*c6 );    
      sum12_3=sum12_3+sum( veh52(i) .*mu2.^(j-1).*mu1.*c6 );



 end

 sum7=v1.*(-sum7_1+sum7_2+sum7_3);
 sum8=b0*(v1^2)*sum8_1+v1*v2*sum8_2+v1^2*sum8_3;
 sum11=b0.*sum11_1+b0.*sum11_2+sum11_3+sum11_4+sum11_5;
 sum12=b0.*v1.*sum12_1+b0.*v1.*2.*sum12_2+v2.*sum12_3;

%compute the inputs appearing in the equation of the third moment
 ST1 = sum(veh3)+3*sum5+3*sum9+6*sum11;       
 ST2 = 3*sum4;                                
 SD1 = sum(veh2)+2*sum1;                      
 ST3 = sum10+2*sum3+2*(sum12+sum7);           
 SD2 = sum(veh1);                             
 ST4 = sum(veh2) + sum1 + sum6 + 2*sum8;      
 SD3 = sum4;                                  

%compute the third moment
 part1 = (T^3) * (r^3);
 part2 = -(3/2)*(T^2)*(r^2)*sum(veh1);
 part3 = 3*T*r*(0.25*SD1+SD2-SD3);
 part4 = -(1/8)*ST1 + ST2 + 0.75*ST3 - (3/2)*ST4;

 er3m=part1+part2+part3+part4;
  
%put results into vectors for output purposes
 vs1m=zeros(13,1);
 vs1m(1)=sum1; vs1m(2)=sum2; vs1m(3)=sum3; vs1m(4)=sum4;   vs1m(5)=sum5;    vs1m(6)=sum6;
 vs1m(7)=sum7; vs1m(8)=sum8; vs1m(9)=sum9; vs1m(10)=sum10; vs1m(11)=sum11;  vs1m(12)=sum12; 
 vs1m(13)=sq5_8;

 vs2m=zeros(7,1);                                 
 vs2m(1)=ST1; vs2m(2)=ST2; vs2m(3)=SD1; vs2m(4)=ST3; vs2m(5)=SD2; vs2m(6)=ST4; vs2m(7)=SD3; 

 vs3m=zeros(4,1);
 vs3m(1)=sum(veh1); vs3m(2)=sum(veh2); vs3m(3)=sum(veh32); vs3m(4)=sum(veh3);




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

