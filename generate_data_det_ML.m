%% setup

%keep it clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%parameters
global I vK vL vCA gK gL gCA vA vB vC vD C phi

I=100; % injected current
vK=-84; % reversal potential for the potassium current (mV)
vL=-60; % reversal potential for the leak current (mV)
vCA=120; % reversal potential for the calcium current (mV)
gK=8; % maximum conductance of the potassium current
gL=2; % leak conductance
gCA=4.4; % maximum conductance of the calcium current
vA=-1.2; % half activation voltage for calcium gating (mv)
vB=18; % reciprocal slope for calcium activation (mv)
vC=2; % half activation voltage for potassium gating (mv)
vD=30; % reciprocal slope for potassium activation (mv)
C=20; % capacitance
phi=0.04; % factor controlling relative timescale of voltage vs n-gate

%time
P=30; % QQQ what is P? (ms)
tmax=85.35*P; % simulation time (ms)
dt=.005; % timestep (ms)


%% system

%initial conditions
x0 = [33.3079; 0.3099569]; % [initial voltage, initial n]

%vector field
F = @(t,u) [1/C*(I-gL*(u(1)-vL)-gK*u(2)*(u(1)-vK)-gCA*m(u(1))*(u(1)-vCA)); ...
    alpha(u(1))*(1-u(2))-beta(u(1))*u(2)];

%solve
[t,U]=ode113(F,0:dt:tmax,x0,opts);

%ease of notation
v = U(:,1); % voltage (mv)
n = U(:,2); % fraction of open potassium gates (dimensionless)

%compute period T
[~,loc]=findpeaks(v,'MinPeakHeight',0);
timerz = t(loc);
T = mean(diff(timerz));

%plot
figure(1)
hold on
plot(v(round(end-T/dt):end),n(round(end-T/dt):end),'k','LineWidth',4);
plot(x0(1),x0(2),'k.','MarkerSize',35)
xlabel('voltage v')
ylabel('n-gate n')
title('Isochrons')
set(gca,'FontSize',15)


%% isochrons

%number of isochrons
num=20;

%store isochrons
count=0;
ONE=nan(num,500);
TWO=nan(num,500);

%phases
phases=0:T/num:T;

%period of the cycle
T=phases(end);

%time step of integration
tau=T/700;

%spatial grid (number of points in isochron)
M=100;

%the number of skipped cycles
k=2;

%forward integration (find the LC)
[~,lc] = ode113(F,0:tau:T,x0,opts);

%spatial resolution (spacing of points in isochron)
dx=(max(lc)-min(lc))'/M;

%center of the limit cycle
center = (max(lc)+min(lc))'/2;

%isochron's initial segment
iso=[x0-M^0.5*dx/10, x0+M^0.5*dx/10];
plot(iso(1,:),iso(2,:),'r.','MarkerSize',30)

%backward integration
for t=0:-tau:-(k+1)*T        
    for i=1:size(iso,2)
        
        %move one step
        iso(:,i)=iso(:,i)-tau*feval(F,t,iso(:,i));
    end
    
    %reset indexing
    i=1;
    
    %remove infinite solutions
    while i<=size(iso,2)
        
        %check if one of the isochron points is too far from
        %the center of the LC
        if any(abs(iso(:,i)-center)>1.5*M*dx)
            
            %if so, remove it
            iso = [iso(:,1:i-1), iso(:,i+1:end)];
        else
            
            %check the other points
            i=i+1;
        end
    end
    
    %reset indexing again...
    i=1;
    
    %adding points on a line ???
    while i<=size(iso,2)-1
        
        %normalized distance
        d=sqrt(sum(((iso(:,i)-iso(:,i+1))./dx).^2));
        
        %add a point in the middle (linear interpolation)
        if d > 2
            iso = [iso(:,1:i), (iso(:,i)+iso(:,i+1))/2 ,iso(:,i+1:end)];
        end
        
        %remove the point
        if d < 0.5
            iso = [iso(:,1:i), iso(:,i+2:end)];
        else
            
            %keep going
            i=i+1;
        end
    end
    
    
    %refresh the screen after each period of integration (clear isochrons)
    if (mod(-t,T)<=tau/2) && (-t<k*T+tau)
        cla;
    end
    
    %plot the isochrons at specified phases
    if min(abs(mod(-t,T)-phases))<tau/2
        
        %update count
        count=count+1;
        
        if (count>1)
            plot(iso(1,:),iso(2,:),'b.');
            drawnow;
            ONE(count,1:length(iso(1,:)))=iso(1,:);
            TWO(count,1:length(iso(1,:)))=iso(2,:);
        end
    end
end

%properly order isochron matrices
ONE=flipud(ONE);
TWO=flipud(TWO);
ONE=ONE(1:num,:);
TWO=TWO(1:num,:);
ONE(ONE == 0) = nan;
TWO(TWO == 0) = nan;

%colormap
mymapboy = jet(num);

%plot the isochrons
for j=1:num
   plot(ONE(j,:),TWO(j,:),'.','color',mymapboy(j,:))
end

%plot
figure(1)
plot(v(round(end-T/dt):end),n(round(end-T/dt):end),'k','LineWidth',4);
set(gca,'FontSize',15)
hold on

%save LC
lc = [v(round(end-T/dt):end)'; n(round(end-T/dt):end)'];


%% write data

%do it
writematrix(lc,'LC.txt');
writematrix(ONE,'isoX.txt');
writematrix(TWO,'isoY.txt');


%% functionals

%alpha
function[u]=alpha(V)
global vC vD phi

xi=(V-vC)/vD;
u=phi*cosh(xi/2)./(1+exp(-2*xi));
end

%beta
function[u]=beta(V)
global vC vD phi

xi=(V-vC)/vD;
u=phi*cosh(xi/2)./(1+exp(2*xi));
end

%m
function[u]=m(V)
global vA vB

xi=(V-vA)/vB;
u=1/2*(1+tanh(xi));
end


%% runtime about 5 sec on a 2023 MacBook Pro (Apple M3 Max)