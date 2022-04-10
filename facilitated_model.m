%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April 2014)

% Facilitated Transport model simulation using optimal set of 
% parameters that were estimated using 'lsqcurvefit' algorithm. 
% Method of lines is used to convert PDEs in the extracellular domain to ODEs



function facilitated_model

% Clear previous files
  clear all
  clc
  
  
% Parameters shared with the ODE routine
   global ncall 
   global wim wme d nc po2 H n c1 c2 c3 c12 kI kM kE lam_im lam_me a v step
   global model rem

   model = input('Please enter a value for model medium alone(1),Cp alone(2),apo-tf alone(3), apo-tf alone+cp(4):');

   if (model==1)
    modeltype='medium alone';
   elseif (model==2)
    modeltype='Cp alone';
   elseif (model==3)
     modeltype='Tf alone';
   else
     modeltype='cp plus tf';
   end
  


  n=41; % Number of horizontal slices in the extracellular domain
  
  nc=11;% total number of species
  
  c1=5; % number of species in the intracellular compartment
  
  c2=10;%number of species in the membrane compartment
  
  c3=7; %number of species in the extracellular domain
  
  c12=c1+c2;
  
  
  
  % volume for EC domain%%%%%%%%%%%%%%% (in ml or cm3)
  v=1;
  
  %initializing mass transfer coefficients and diffusivities, 11 species
  
  d=zeros(11,1);
  lam_im=zeros(11,1);
  lam_me=zeros(11,1);
  

  %%%% Optimal rate constant values 
  
   k1 = 32.5; % (per (M.s) %ferrous oxidation fwd by o2
   k11= 121; % (per s)    %ferric reduction
   k2 = 88; % (per (M.s)  %ferrous fpn formation
   k3 = 3.86; % (per (M.s)  % ferrous fpn oxidation by cp (CP ACTION)
   k4 = 15.397; % (per (M.s)  % TF ACTION (EXTRACELLULAR)
   k5 = 18.6;%   (per (M.s)) % CP ACTION (EXTRACELLULAR)
   k6 = 3.4; % (per s) (from manuscript)% fwd disscoiation of ferric-fpn
   k66= 1.7281; % (per (M.s))% reversible formation of ferric from ferric-fpn
   rem=0.0818; %rate of sequestering ferrous to ferritin %kr in the text
   
  
    
  kI=[0;0;k2];
  kM=[k1;k11;k3;k6;k66;k4];
  kE=[k5;k4];
  

  %specify diffusivities%(in cm^2/sec)
  dfe=3.46*10^-6;
  dtf=9.02*10^-8;
  dcp=7.13*10^-9;
  do2=3.24*10^-6;
 
  d(1)=dfe;
  d(2)=dfe;
  d(3)=0.1*dfe; % diffusivity of FPN (and also associated species) is 1/10th of that of ferrous ion
  d(4)=0.1*dfe;
  d(5)=do2;
  d(6)=0.1*dfe;
  
  fe0=2.47*10^-6;% in M
  tf0=55*10^-6; % in M
  cp0=3*10^-6; % in M
 
  % all values in per second 
  wim=[0.297;0.01;0.00391];%mass transfer coefficient from i to m
  wme=[0.09;0.00391;0.0329;0.0329;0.277;0.277;0.277];%mass transfer coefficient from m to e
  
  %mass transfer coefficients between compartments

 lam_im(3)=wim(1);
 lam_im(4)=wim(2);
 lam_im(5)=wim(3);

 lam_me(2)=wme(1);
 lam_me(5)=wme(2);
 lam_me(7)=wme(3);
 lam_me(8)=wme(4);
 lam_me(9)=wme(5);
 lam_me(10)=wme(6);
 lam_me(11)=wme(7);
  
  
  po2=0.01;%partial pressure of oxygen 
  H=769.23*10^-6;%henry's constant (atm/uM)
  
  zl=0.0;%lower limit of z
  zu=1.0;%upper limit of z
  dz=(zu-zl)/(n-1);
  
  
  switch (model)
      case(1)
          kM(3)=0;% no Cp
          kM(6)=0;% no Tf
          kE(1)=0;% no cp
          kE(2)=0;% no tf
      case(2)
          kM(6)=0;% no Tf
          kE(2)=0;% no Tf
          d(7)=dcp;
          d(8)=dcp;
      case(3)
          kM(3)=0; % no Cp
          kE(1)=0; % no Cp
          d(9)=dtf;
          d(10)=dtf;
          d(11)=dtf;
      case(4)
          d(7)=dcp;
          d(8)=dcp;
          d(9)=dtf;
          d(10)=dtf;
          d(11)=dtf;
  end 
  
  
% Initial condition all concentration stored in one-dimensional u matrix  
  
% TOTAL 11 SPECIES , (5 intra, 10 mem, 7 extra )

 u=zeros(c12+c3*n,1);%   
 
 %initialize concentrations in the intracellular domain
 u(1,1)=fe0;%ferrous concentration inside the cell
 u(3,1)=fe0*2;%FPN concentration inside the cell
 u(5,1)=po2/H;%oxygen concentration inside the cell
 
 %initialize concentrations in the membrane domain
 u(9,1)=po2/H;%oxygen concentration inside the membrane domain, rest all species are zero at start 
 
 
 
%Specify the initial conditions depending on experimental condition selected
 
 for i=1:n
     u(c12+(i-1)*c3+2,1)=po2/H;
     %u(c12+(i-1)*c3+1,1)=fe0/100;%ferric concentration outside the cell
     if (model==2)%only cp added
        u(c12+(i-1)*c3+4,1)=cp0;%cu+2 
     elseif (model==3)%only tf added
        u(c12+(i-1)*c3+5,1)=tf0;%tf 
     elseif (model==4) % both cp and tf added
        u(c12+(i-1)*c3+4,1)=cp0;%cu+2 
        u(c12+(i-1)*c3+5,1)=tf0;%tf 
     end
 
 end
 
 
 


% Independent variable for ODE integration
  t0=0.0;
  tf=60*50;
  step=1;
  tout=t0:1:tf; 
  nout=tf/step;
  ncall=0;
%
% ODE itegration
  reltol=1.0e-02; abstol=1.0e-02;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
   
  
  u0=u(:,1);
  length(u0');
  
  
  [t,u]=ode15s(@system,tout,u0,options);
  
  
  fet=zeros(tf/step,1);
  fet_ferric=zeros(tf/step,1);
  per_fet_ferric=zeros(tf/step,1);
  
  % Calculating total iron released in the extracellular medium (that is 'fet") using simpsons rule of integration  

 
for j=1:tf/step
    
        for i=1:n
            if (i==1)
                fet(j)=fet(j)+u(j,c12+(i-1)*c3+1)+u(j,c12+(i-1)*c3+6)+2*u(j,c12+(i-1)*c3+7);
                fet_ferric(j)=fet_ferric(j)+u(j,c12+(i-1)*c3+1);
                fet;
                
            elseif (i==n)
                
                fet(j)=fet(j)+u(j,c12+(i-1)*c3+1)+u(j,c12+(i-1)*c3+6)+2*u(j,c12+(i-1)*c3+7);
                fet_ferric(j)=fet_ferric(j)+u(j,c12+(i-1)*c3+1);
                fet;  
            elseif (mod(i-1,2)==1) % even points, multiplied by 2
                
                fet(j)=fet(j)+2*(u(j,c12+(i-1)*c3+1)+u(j,c12+(i-1)*c3+6)+2*u(j,c12+(i-1)*c3+7));
                fet_ferric(j)=fet_ferric(j)+2*u(j,c12+(i-1)*c3+1);
                fet;
            else
                % odd points, multiplied by 4
                fet(j)=fet(j)+4*(u(j,c12+(i-1)*c3+1)+u(j,c12+(i-1)*c3+6)+2*u(j,c12+(i-1)*c3+7));
                fet_ferric(j)=fet_ferric(j)+4*u(j,c12+(i-1)*c3+1);
                fet;
            end
        end
        
        fet;
        
      
      
      fet(j)=fet(j)*dz/3; % applying final simpsons rule (divide by 3) to get iron in molar
      
        
      
      fet_ferric(j)=fet_ferric(j)*dz/3;
      
      fet(j)=fet(j)*v/1000;% getting iron concentration in moles
        
      fet_ferric(j)=fet_ferric(j)*v/1000;
      
      per_fet_ferric(j)=fet_ferric(j)/fet(j)*100;
      
        %error(j)=(abs(fet(j)*1000/v-(fe0-u(j,1))))*100/(fe0-u(j,1));
%         fprintf('\n time = %2d   fet = %2d fet_ferric = %2d fe0-u(j,1) = %2d error = %2d gradient = %5d\n ',...
%         tout(j)/60, fet(j), fet_ferric(j), fe0-u(j,1), error(j), (u(j,7)-u(j,4)));
%     
    if (tout(j)/60==30) % iron release at 30 min
        fet(j)
    end
   
end
  


figure(3)% bar graph of iron release
bar(model,fet(tf/step),'r');
hold on;

%dlmwrite('fet_apocp.new2.csv',fet);


figure(5)% plot the total iron release
plot(fet(:)*10^12, '-k')
hold on;



%%function for defining odes where dss044 is called

function ut=system(t,u)

global ncall n nc  d  c1 c2 c12 c3 lam_im lam_me 



% defining the ode/PDE system
%dividing u values into ui, um and ue for the I, M and E compartments


ui=zeros(nc,1);
um=zeros(nc,1);
ue=zeros(n,nc);

for i=1:c1
    ui(i)=u(i);
    
end



for i=1:c2
    um(i+1)=u(i+c1);
end 



for i=1:n
    ue(i,2)=u((i-1)*c3+c12+1);
    ue(i,5)=u((i-1)*c3+c12+2);
    ue(i,7)=u((i-1)*c3+c12+3);
    ue(i,8)=u((i-1)*c3+c12+4);
    ue(i,9)=u((i-1)*c3+c12+5);
    ue(i,10)=u((i-1)*c3+c12+6);
    ue(i,11)=u((i-1)*c3+c12+7);
end



k=1;
  
  if (le(k,(c12+c3*n)))
   dummy=0;
   dummy1=0;
   dummy2=0;
        for i=1:c1
%        intracellular compartment odes
                
                    dummy=i;
                    dummy1=k;
                    ri=reacti(ui);
                    ut(dummy1,1)=-lam_im(dummy)*(ui(dummy)-um(dummy))+ri(dummy);
                    if (i==3)
                        if (ui(3)>um(3))
                            ut(dummy1,1)=ri(dummy);
                        end 
                    elseif (i==4)
                        if (ui(4)<um(4))
                            ut(dummy1,1)=ri(dummy);
                        end
                    end
 
                    k=k+1;
                    
        end 
        
       
        dummy=0;
        dummy1=0;
        dummy2=0;
        
        
        for i=2:nc
%        membrane compartment odes
                
                    dummy=i;
                    dummy1=k;
                    rm=reactm(um);
                    ut(dummy1,1)=lam_im(dummy)*(ui(dummy)-um(dummy))-lam_me(dummy)*(um(dummy)-ue(1,dummy))+rm(dummy);
                   
                    if (and((i==3),(ui(3)>um(3))))
                            ut(dummy1,1)=rm(dummy); 
                    elseif (and((i==4),(ui(4)<um(4))))
                            ut(dummy1,1)=rm(dummy);

                    end
 
                    k=k+1;
                    
        end 
        
        
        dummy=0;
        dummy1=0;
        dummy2=0;
       % extracellular domain final odes (pdes are converted to odes using
       % method of lines , using dss routine
        
       for i=1:n
                dummy=i;  
               
                
            for j=1:nc
                dummy=i; %space
                dummy1=k; % time
                dummy2=j; % species
                re=reacte(ue);
                uxx=dss(ue,um);
                
                    
                switch(j)
                    
                    case(2)
                        %%%convert diffusivty in cm2 to mm2 and dz is 1/40th of 1mm
                        
                        ut(dummy1,1)=(d(dummy2)/(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                    case(5)
                        
                        ut(dummy1,1)=(d(dummy2)/(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                    case(7)
                        
                        ut(dummy1,1)=(d(dummy2)/(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                    case(8)
                        
                        ut(dummy1,1)=(d(dummy2)/(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2); 
                        
                        
                        k=k+1;
                     case(9)
                        
                        ut(dummy1,1)=(d(dummy2)*(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                     case(10)
                        
                        ut(dummy1,1)=(d(dummy2)*(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                     case(11)
                        
                        ut(dummy1,1)=(d(dummy2)*(0.25^2))*uxx(dummy,dummy2)+re(dummy,dummy2);
                        
                        
                        k=k+1;
                      
                
                end
        
               
            end  
            
        end
          
        
  else
      
      
  end

 
 % Increment calls to this function
   ncall=ncall+1;
   
   
 
 length(ut);
    
   
%defining the intracellular reaction system%%%%%%%%%

function ri=reacti(ui)
%%

global  rem kI 


ri(1)=-kI(1)*ui(1)*ui(5)-kI(3)*ui(1)*ui(3)+kI(2)*ui(2)-rem*ui(1);
ri(2)=kI(1)*ui(1)*ui(5)-kI(2)*ui(2);
ri(3)=-kI(3)*ui(1)*ui(3);
ri(4)=kI(3)*ui(1)*ui(3);
ri(5)=-kI(1)*ui(1)*ui(5)+kI(2)*ui(2);

  



%defining the membrane reaction system%%%%%%%%%

function rm=reactm(um)

global   kM 



rm(1)=0;%ferrous
rm(2)=kM(4)*um(6)-kM(5)*um(2)*um(3)-kM(6)*um(2)*um(9)-kM(6)*um(2)*um(10);%ferric
rm(3)=kM(4)*um(6)-kM(5)*um(2)*um(3);%FPN
rm(4)=-kM(1)*um(4)*um(5)+kM(2)*um(6)-kM(3)*um(4)*um(8);%ferrous-fpn
rm(5)=(-kM(1)*um(4)*um(5)+kM(2)*um(6))/2;%o2
rm(6)=kM(1)*um(4)*um(5)-kM(2)*um(6)+kM(3)*um(4)*um(8)-kM(4)*um(6)+kM(5)*um(2)*um(3);%ferric-fpn
rm(7)=kM(3)*um(4)*um(8);%cu+1
rm(8)=-kM(3)*um(4)*um(8);%cu+2
rm(9)=-kM(6)*um(2)*um(9);%Tf
rm(10)=kM(6)*um(2)*um(9)-kM(6)*um(2)*um(10);%monoferric-tf
rm(11)=kM(6)*um(2)*um(10);%diferric-tf
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining the extracellular reaction system%%%%%%%%%
function re=reacte(ue)

global  n  nc   kE 

re=zeros(n,nc);
    
for i=1:n
re(i,2)=-kE(2)*ue(i,2)*ue(i,9)-kE(2)*ue(i,2)*ue(i,10);%ferric
re(i,5)=-kE(1)*ue(i,5)*ue(i,7);%o2
re(i,7)=-kE(1)*ue(i,5)*ue(i,7);%cu+1
re(i,8)=kE(1)*ue(i,5)*ue(i,7);%cu+2
re(i,9)=-kE(2)*ue(i,2)*ue(i,9);%Tf
re(i,10)=kE(2)*ue(i,2)*ue(i,9)-kE(2)*ue(i,2)*ue(i,10);%monoferric-tf
re(i,11)=kE(2)*ue(i,2)*ue(i,10);%diferric-tf
end
  

  
%%Application of the method of lines to approximate the second derivatives in the PDE for extracellular domain,%thereby reducing the PDE system to ODE system
function uxx=dss(ue,um)

% Problem parameters
  global n model d nc po2 H   lam_me
  
  uxx=zeros(n,nc);
   
      xl=0.0;
      xu=1.0;
      
 for i=1:nc
     
      u1=zeros(n,1);
      u1x=zeros(n,1);


    dumb=i;
      for j=1:n
        u1(j,1)=ue(j,dumb);
      end
      
      nl=2; % Neumann
      nu=2; % Neumann
    
 % BC at x = 0   
      switch(dumb)
          
          %%%"lam_me" are "beta_me" as per manuscript 
        case(2)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);
        case(5)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);   
        case(7)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);
        case(8)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);   
        case(9)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);
        case(10)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);
        case(11)     
           u1x(1)=-lam_me(dumb)*(um(dumb)-u1(1))/d(dumb);   
      end
      
% BC at x = 1
      if (dumb==5)
        u1x(n)=lam_me(5)*((po2/H)-u1(n))/d(5); 
      end
  
     
% Calculate uxx using second derivative routines
  
  u1xx=dss044(xl,xu,n,u1,u1x,nl,nu); % second order
  
  u1xx=u1xx';
  for z=1:n
  uxx(z,dumb)=u1xx(z,1);
  end


  
  switch(model)
      case(1)
          if(i>5)
              uxx(:,i)=0;
          end 
      case(2)
          if (or(i==9,i==10))
              uxx(:,i)=0;
          elseif(i==11)
              uxx(:,i)=0;
          end
      case(3)
          if(or(i==7,i==8))
              uxx(:,i)=0;
          end
  end
 end
  
 
              
