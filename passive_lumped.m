%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April 2014)

% Passive gradeint (spatially lumped version) model simulation using optimal set of 
% parameters that were estimated using 'lsqcurvefit' algorithm
 
function passive_lumped

% Clear previous files
  clear all
  clc
  
  
% Parameters shared with the ODE routine

   global ncall 
   global w po2 H tf0 cp0
   global model 

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
  


%Initial concentration of species

  fe0=2.47*10^-6;% in M
  
  tf0=55*10^-6;% in M
  
  cp0=3*10^-6;% in M. 
  
  w=3.75*10^-5;%mass transfer coefficient for ferrous
  
  po2=0.01;%partial pressure of oxygen
 
  H=769.23;%henry's constant (atm/M)
  
 
  u=zeros(7,1);% all concentration stored in one-dimensional u matrix 
 
  u(1,1)=fe0;
 
  u(7,1)=po2/H;
  
  
 
% Independent variable for ODE integration
  t0=0.0;
  
  tf=60*30; %Iron release time interval in s
  
  step=1;
  
  tout=linspace(t0,tf,tf/step); 
  
  nout=tf/step;
  
  ncall=0;
%
% ODE itegration
  reltol=1.0e-04; abstol=1.0e-04;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
   
  
  u0=[fe0;0;0;0;0;0;po2/H]; %initial concentration of species
  
  
  [t,u]=ode15s(@extra,tout,u0,options);
  
  fet=zeros(tf/step,1);
  
  error=zeros(tf/step,1);
  
 %calculating total iron release in the extracellular medium ('fet") (in molar) 
 
        for j=1:tf/step
            fet(j)=fet(j)+u(j,2)+u(j,3)+u(j,4)+2*u(j,5);
            fet;
        end
       
        fet(j)=fet(j)*1/1000;% iron release  in moles
        
        error(j)=(abs(fet(j)-(fe0-u(j,1))))*100/(fe0-u(j,1));
        fprintf('\n time = %2d   fet = %2d fe0-u(j,1) = %2d error = %2d gradient = %5d\n ',...
        tout(j), fet(j), fe0-u(j,1), error(j), w*(u(j,2)-u(j,1)));

  
  figure(2)
  
  plot(fet(tf/step),fe0-u(tf/step,1),'go');

  figure(3)% bar graph of iron release

  bar(model,fet(tf/step),'r');

  hold on;


  figure(5)% plot the total iron release

  plot(fet(:)*10^12, 'b-')

  hold on;



function ut=extra(t,u)

global ncall w model po2 H tf0 cp0

% defining the ode system

% Parameters 

        k1 = 32.5; % (per (M.s)
        k11= 121; % (per s)
        k2 = 3.86; % (per (M.s)
        k3 = 18.6; % (per (M.s)
        k4 = 15.4; % (per (M.s)
        kr = 0.08; % (per s) 

    cp=cp0-u(7); % Equation for Cu2+ ions
    
    trf=tf0-u(4)-u(5); % Equation for transferrin
    
  
%Specify the initial conditions depending on experimental condition selected    
    if (model==1)
        k2=0.0;
        k3=0.0;
        k4=0.0;
    elseif (model==2)
        k4=0.0;
    elseif (model==3)
        k2=0.0;
        k3=0.0;
    elseif (model==4)
        k1 = 32.5; % (per (M.s)
        k11= 121; % (per s)
        k2 = 3.86; % (per (M.s)
        k3 = 18.6; % (per (M.s)
        k4 = 15.4; % (per (M.s)
    end
    
    
  ut=zeros(7,1);
  
  
  
  ut(1)=w*(u(2)-u(1))-kr*u(1,1); % IC Fe2+ 
  
  ut(2)=-w*(u(2)-u(1))-k1*u(2)*u(7)+k11*u(3)-k3*u(2)*cp; % EC Fe2+
  
  ut(3)=k1*u(2)*u(7)-k11*u(3)+k3*u(2)*cp-k4*u(3)*trf-k4*u(3)*u(4); % Fe3+
  
  ut(4)=k4*u(3)*trf-k4*u(3)*u(4);% Monoferric Tf
  
  ut(5)=k4*u(3)*u(4); % Holo-Tf
  
  ut(6)=-k2*u(6)*u(7)+k3*u(2)*cp; % Cp1+
 
  ut(7)=-k1*u(2)*u(7)+k11*u(3)-k2*u(6)*u(7)+w*((po2/H)-u(7)); %Oxygen
  
 

  
 % Increment calls to this function
   ncall=ncall+1;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
