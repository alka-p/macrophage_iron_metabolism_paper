%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April
% 2014)

% Passive gradient (distributed version) model simulation using optimal set of 
% parameters that were estimated using 'lsqcurvefit' algorithm. 
% Method of lines is used to convert PDEs in the extracellular domain to ODEs


function passive_model_distributed


% Clear previous files
   clear all
   clc
  
  
% Parameters shared with the ODE routine
   global ncall 
   global w d nc po2 H n
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


  n=41; % Number of horizontal slices in the extracellular domain 
        
  nc=8; % Number of pde species 
  
  
  %specify diffusivities %(in mm^2/sec)
  
  dfe=3.46*10^-4; % Iron species diffusivity
  
  dtf=9.02*10^-6; % Tf related species diffusivity
  
  dcp=7.13*10^-6; % Cp related species diffusivity
  
  do2=3.24*10^-4; % O2 diffusivity
  
  d(1)=dfe;
  d(2)=dfe;
  d(3)=dtf;
  d(4)=dtf;
  d(5)=dtf;
  d(6)=dcp;
  d(7)=dcp;
  d(8)=do2;
  
  % volume for EC domain%%%%%%%%%%%%%%% (in ml or cm3)
  
  v=1;

  %Initial concentration of species

  fe0=2.47*10^-6;%in M

  tf0=55*10^-6;%in M
  
  cp0=3*10^-6; %in M. 
  
  w=3.75*10^-5;%mass transfer coefficient for ferrous
  
  po2=0.01;%partial pressure of oxygen

  H=769.23;%henry's constant (atm/M)
  
  zl=0.0;%lower limit of z
  zu=1.0;%upper limit of z in mm
  dz=(zu-zl)/(n-1);
  
  

% All concentrations stored in ond-dimensional u matrix  
% TOTAL NINE SPECIES (ODE +PDE),[fe2+(IC) fe+2(EC) fe+3 mono-tf holo-tf tf cp+1 cp+2 o2 ];

 u=zeros(1+nc*n,1);% example for 20 grids , 153 odes (1+nc*(n-1))  
 
 u(1,1)=fe0;
 
 

 
 for i=1:n
     u(i*nc+1,1)=po2/H;
 end
 
 %Specify the initial conditions depending on experimental condition selected
 
 if (model==2)%only cp added
      for i=0:n-1
          u(i*nc+8,1)=cp0;%cu+2 is at 8,16,24...
      end 
 elseif (model==3)%only tf added
      for i=0:n-1
          u(i*nc+6,1)=tf0;%tf at 6,14,22....
      end
 elseif (model==4)     
      for i=0:n-1
         u(i*nc+8,1)=cp0;
         u(i*nc+6,1)=tf0;
      end
 end
 
 
 

% Independent variable for ODE integration
  t0=0.0;
  tf=60*30;
  step=1;
  tout=linspace(t0,tf,tf/step); 
  nout=tf/step;
  ncall=0;

% ODE itegration
  reltol=1.0e-04; abstol=1.0e-04;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
   
  
  u0=u(:,1); %length(u0')
  
  [t,u]=ode15s(@extra,tout,u0,options);
  
  fet=zeros(tf/step,1);
  error=zeros(tf/step,1);
  
 
  
% Calculating total iron released in the extracellular medium (that is 'fet") using simpsons rule of integration  

for j=1:tf/step
    
        for i=0:n-1
            if (i==0)
                fet(j)=fet(j)+u(j,2+nc*i)+u(j,3+nc*i)+u(j,4+nc*i)+2*u(j,5+nc*i);
                fet;
                
            elseif (i==n-1)
                
                fet(j)=fet(j)+u(j,2+nc*i)+u(j,3+nc*i)+u(j,4+nc*i)+2*u(j,5+nc*i);
                fet;   
            elseif (mod(i,2)==1)
                
                fet(j)=fet(j)+2*(u(j,2+nc*i)+u(j,3+nc*i)+u(j,4+nc*i)+2*u(j,5+nc*i));
                fet;
            else
                
                fet(j)=fet(j)+4*(u(j,2+nc*i)+u(j,3+nc*i)+u(j,4+nc*i)+2*u(j,5+nc*i));
                fet;
            end
        end
        
        
        fet(j)=fet(j)*dz/3;
        fet(j)=fet(j)*v/1000;%  fet  in moles
        
        error(j)=(abs(fet(j)-(fe0-u(j,1))))*100/(fe0-u(j,1));
        fprintf('\n time = %2d   fet = %2d fe0-u(j,1) = %2d error = %2d gradient = %5d\n ',...
        tout(j), fet(j), fe0-u(j,1), error(j), abs(w*(u(j,2)-u(j,1))));
end
  
figure(3)% bar graph of iron release
bar(model,fet(tf/step),'r');
hold on;



figure(5)% plot the total iron release in pmoles
plot(fet(:)*10^12, 'b-')
hold on;


%%function for defining odes where dss044 is called

function ut=extra(t,u)

global ncall d w nc n

% defining the ode/PDE system
   
  kr=0.0818; % (per s)
  
  ut(1,1)=w*(u(2,1)-u(1,1))-kr*u(1,1); % ODE for intracellular domain


  k=2;
  
  if (le(k,(1+nc*n)))
   
        for i=1:n
%          for j=1:nc
%             u1(j)=u(1+(i-1)*nc+j); 
%          end
                for j=1:nc
                    dummy=j;
                    dummy1=k;
                    dummy2=i;
                    r=react(dummy2,u);
                    uxx_new=dss(dummy,u);
                    uxx_new(dummy2);
                    ut(dummy1,1)=(d(dummy)/(0.25^2))*uxx_new(dummy2)+r(dummy);
                    k=k+1;
                    length(ut);
                end 
        end
  else
      
  end

  
 % Increment calls to this function
   ncall=ncall+1;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining the reaction system%%%%%%%%%

function r=react(ni,u1)

global model nc n

%% defining the parameters 

        k1 = 32.5; % (per (M.s)
        k11= 121; % (per s)
        k2 = 3.86; % (per (M.s)
        k3 = 18.6; % (per (M.s)
        k4 = 15.4; % (per (M.s)

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
    
        u=zeros(nc*n,1);
  
        for j=1:nc
            u(j)=u1(1+(ni-1)*nc+j);
        end   
    
       
    r(1)=-k1*u(1)*u(8)+k11*u(2)-k3*u(1)*u(7);
    
    r(2)=k1*u(1)*u(8)-k11*u(2)+k3*u(1)*u(7)-k4*u(2)*u(5)-k4*u(2)*u(3);
    
    r(3)=k4*u(2)*u(5)-k4*u(2)*u(3);

    r(4)=k4*u(2)*u(3);
 
    r(5)=-k4*u(2)*u(5);

    r(6)=-k2*u(6)*u(8)+k3*u(1)*u(7);
 
    r(7)=k2*u(6)*u(8)-k3*u(1)*u(7);
 
    r(8)=-k1*u(1)*u(8)+k11*u(2)-k2*u(6)*u(8);
  
   r=r';
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
%%Application of the method of lines to approximate the second derivatives in the PDE,thereby reducing the PDE system to ODE system
function uxx_new=dss(nj,u)


  global d w nc n po2 H
  
  
  %defining new vector of concentrations without the the intracellular concentration, that is only extracellular concentrations
  
      xl=0.0;
      xu=1.0;
      u2=zeros(n,1);
      u2x=zeros(n,1);
      
      for j=0:n-1 
         u2(j+1,1)=u(j*nc+(nj+1));
      end
      
      
       
% BC at x = 0
      if (nj==1)
        u2x(1)=w*(u(2,1)-u(1,1))/d(1);
      end
      
% BC at x = 1
      if (nj==8)
        u2x(n)=w*((po2/H)-u(1+nc*n,1))/d(8);
      end
  
 
% Calculate uxx using second derivative routines

  nl=2; % Neumann
  nu=2; % Neumann
  
  u2xx=dss044(xl,xu,n,u2,u2x,nl,nu); % second order
  
  uxx_new=u2xx;
  
