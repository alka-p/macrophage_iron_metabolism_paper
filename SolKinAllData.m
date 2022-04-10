%Author: Alka Potdar
% Run using MATLAB R2009b
% Computational modeling and analysis of macrophage iron release (April 2014)

% Fitting to entire solution kinetics experimental data that is compute k1,k11 k4,k6,k2,k3
% simulataneously using both noCp and withCp data sets

function SolKinAllData
clear 
clc

 
% Experimental data without cp

xdata1=importdata('xdata_nocp.mat');
ydata1=importdata('ydata_nocp.mat');

% Experimental data without cp

xdata2=importdata('x_new_withcp.mat');
ydata2=importdata('y_new_withcp.mat');

% Merging the two datasets 

xdata=[xdata1;xdata2];
ydata=[ydata1;ydata2];


% ydata is in millimolar, convert to micromolar

xdata=transpose(xdata);
ydata=ydata*1000;


% Initial values of k1,k11,k4,k6,k2 and k3 found using individual datasets  followed by multiple lsqcurvefit runs

x0=[0.0044  163.2387    0.0095    1.3430    2.9163    8.4519];

% 
% %combined least squares fit
[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@func,x0,xdata,ydata,[0.00001 0.00001 0.00001 0.00001 0.00001 0.00001],[1000 1000 1*1000 1*1000 1*1000 1*1000])



%--------------------------------------------------------------------------

% function, f for iron transferrin (ferric bound to transferrin, monoferric and holoferric forms)

%--------------------------------------------------------------------------
% computing the values of function f using ode15s


function f=func(x,xdata)

xdata1=xdata(1,1:25);

xdata2=xdata(1,26:74);

y1 = ode15s(@odefun1,[0 120],[120 0 0 0 55],[],x(1),x(2),x(3),x(4));

y2 = ode15s(@odefun2,[0 120],[120 0 0 0 55 0 0.9],[],x(1),x(2),x(3),x(4),x(5),x(6));

% Fe+2: 120 uM; Tf: 55 uM; Cp2+ : 0.9 uM 

sxint1 = deval(y1,xdata1);

f1=(sxint1(3,:)+2*sxint1(4,:));

sxint2 = deval(y2,xdata2);

f2=(sxint2(3,:)+2*sxint2(4,:));


f=[transpose(f1);transpose(f2)];

dlmwrite('f.csv',f);

dlmwrite('xsim_nocp_both.csv',xdata1);

dlmwrite('ysim_nocp_both.csv',f1);


dlmwrite('ysim_withcp_both.csv',f2);

dlmwrite('xsim_withcp_both.csv',xdata2);

return


%--------------------------------------------------------------------------
% defining the set of ordinary differential equations for the reaction  system without Cp

function df1=odefun1(t,c,k1,k11,k4,k6)


df1 = zeros(5,1);    % a column vector

 
co2= 220; %concentration of oxygen (uM) for 1% oxygen
    
df1(1)=-k1*c(1)*co2+k11*c(2); % eqn for fe+2

df1(2)=k1*c(1)*co2-k11*c(2)-k4*c(2)*c(5)-k4*c(2)*c(3)-k6*c(2); % eqn for fe+3

df1(3)=k4*c(2)*c(5)-k4*c(2)*c(3); % eqn for monoferric transferrin

df1(4)=k4*c(2)*c(3); % eqn for holo transferrin

df1(5)=-k4*c(2)*c(5); % eqn for transferrin

return




%--------------------------------------------------------------------------
% defining the set of ordinary differential equations for the reaction  system with Cp
function df2=odefun2(t,c,k1,k11,k4,k6,k2,k3)

df2 = zeros(7,1);    % a column vector

co2= 220; %concentration of oxygen (uM) for 1% oxygen
    

df2(1)=-k1*c(1)*co2+k11*c(2)-k3*c(1)*c(7);% eqn for fe+2

df2(2)=k1*c(1)*co2-k11*c(2)+k3*c(1)*c(7)-k4*c(2)*c(5)-k4*c(2)*c(3)-k6*c(2);% eqn for fe+3,k6 isfirst order decay rate constant for ferric ions in solution

df2(3)=k4*c(2)*c(5)-k4*c(2)*c(3);% eqn for monoferric transferrin

df2(4)=k4*c(2)*c(3);% eqn for holo transferrin

df2(5)=-k4*c(2)*c(5);% eqn for transferrin

df2(6)=-k2*c(6)*co2+k3*c(1)*c(7);%eqn for cp+1

df2(7)=k2*c(6)*co2-k3*c(1)*c(7);%eqn for cp+2




%--------------------------------------------------------------------------


return























