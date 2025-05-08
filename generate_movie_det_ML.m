%% setup

%mr clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%parameters
global I vK vL vCA gK gL gCA vA vB vC vD C phi

I=100;
vK=-84;
vL=-60;
vCA=120;
gK=8;
gL=2;
gCA=4.4;
vA=-1.2;
vB=18;
vC=2;
vD=30;
C=20;
phi=0.04;

%vector field
F = @(t,u) [1/C*(I-gL*(u(1)-vL)-gK*u(2)*(u(1)-vK)-gCA*m(u(1))*(u(1)-vCA)); ...
    alpha(u(1))*(1-u(2))-beta(u(1))*u(2)];

%time
P=1;
tmax=85.35*P;
dt=.005;
t = 0:dt:tmax;

%number of isochrons
num = 20;


%% load in data and make initial plot

%LC
LC = load('LC.txt');
isoX = load('IsoX.txt');
isoY = load('IsoY.txt');

%colormap
mymapboy = jet(num);

%plot the isochrons
for j=1:num
    figure(1)
    plot(isoX(j,:),isoY(j,:),'.','color',mymapboy(j,:))
    hold on
end

%plot
figure(1)
plot(LC(1,:),LC(2,:),'k','LineWidth',4);
set(gca,'FontSize',15)
axis square
box on
xlim([-100 100])
ylim([0 0.9])


%% prune the isochrons and make plot of initial conditions

%prune
isoXX = isoX(:,110:5:205);
isoYX = isoY(:,110:5:205);

%plot the pruned isochrons
for j=1:num
    figure(1)
    plot(isoXX(j,:),isoYX(j,:),'.','color',mymapboy(j,:),'markersize',20)
    hold on
end


%% solve

%make solution vector
SOL = zeros(num*2*20,length(t));

%count
count = 1;

%loop over total number of isochrons
for i=1:num
    i
    %loop over each point on the isochron
    for j = 1:20

        %find initial condition
        x0 = [isoXX(i,j);isoYX(i,j)];

        %solve
        [~,U]=ode113(F,t,x0,opts);    %%%%% perhaps using the euler-maruyama here with a small amount of noise would illustrate things for the stochastic case?

        %insert into SOL
        SOL(count:count+1,:) = U';

        %update count
        count = count+2;
    end
end

%save the solution trajectories
writematrix(SOL,'solutions.txt');


%% make the movie

%load in data
SOL = load('solutions.txt');

%cheat
onex = SOL(1:2:40,:);
oney = SOL(2:2:40,:);
twox = SOL(41:2:80,:);
twoy = SOL(42:2:80,:);
threex = SOL(81:2:120,:);
threey = SOL(82:2:120,:);
fourx = SOL(121:2:160,:);
foury = SOL(122:2:160,:);
fivex = SOL(161:2:200,:);
fivey = SOL(162:2:200,:);
sixx = SOL(201:2:240,:);
sixy = SOL(202:2:240,:);
sevenx = SOL(241:2:280,:);
seveny = SOL(242:2:280,:);
eightx = SOL(281:2:320,:);
eighty = SOL(282:2:320,:);
ninex = SOL(321:2:360,:);
niney = SOL(322:2:360,:);
tenx = SOL(361:2:400,:);
teny = SOL(362:2:400,:);

onexx = SOL(401:2:440,:);
oneyx = SOL(402:2:440,:);
twoxx = SOL(441:2:480,:);
twoyx = SOL(442:2:480,:);
threexx = SOL(481:2:520,:);
threeyx = SOL(482:2:520,:);
fourxx = SOL(521:2:560,:);
fouryx = SOL(522:2:560,:);
fivexx = SOL(561:2:600,:);
fiveyx = SOL(562:2:600,:);
sixxx = SOL(601:2:640,:);
sixyx = SOL(602:2:640,:);
sevenxx = SOL(641:2:680,:);
sevenyx = SOL(642:2:680,:);
eightxx = SOL(681:2:720,:);
eightyx = SOL(682:2:720,:);
ninexx = SOL(721:2:760,:);
nineyx = SOL(722:2:760,:);
tenxx = SOL(761:2:800,:);
tenyx = SOL(762:2:800,:);


%% finally do it

%make the object
v = VideoWriter('isochrons.mp4','MPEG-4');
open(v)

%loop over each solution time
for k=1:50:length(SOL)
        
    %plot LC
    figure(1)
    plot(LC(1,:),LC(2,:),'k','LineWidth',4);
    set(gca,'FontSize',15)
    axis square
    box on
    xlim([-100 100])
    ylim([0 0.9])
    hold on
    xlabel('v')
    ylabel('n')
    title('Morris Lecar System')
    
    %plot the isochrons
    for j=1:num
        figure(1)
        plot(isoX(j,:),isoY(j,:),'.','color',mymapboy(j,:))
        hold on
    end
    
    %plot isochrons things
    figure(1)
    hold on
    plot(onex(:,k),oney(:,k),'.','color',mymapboy(1,:),'markersize',20)
    plot(twox(:,k),twoy(:,k),'.','color',mymapboy(2,:),'markersize',20)
    plot(threex(:,k),threey(:,k),'.','color',mymapboy(3,:),'markersize',20)
    plot(fourx(:,k),foury(:,k),'.','color',mymapboy(4,:),'markersize',20)
    plot(fivex(:,k),fivey(:,k),'.','color',mymapboy(5,:),'markersize',20)
    plot(sixx(:,k),sixy(:,k),'.','color',mymapboy(6,:),'markersize',20)
    plot(sevenx(:,k),seveny(:,k),'.','color',mymapboy(7,:),'markersize',20)
    plot(eightx(:,k),eighty(:,k),'.','color',mymapboy(8,:),'markersize',20)
    plot(ninex(:,k),niney(:,k),'.','color',mymapboy(9,:),'markersize',20)
    plot(tenx(:,k),teny(:,k),'.','color',mymapboy(10,:),'markersize',20)
    
    plot(onexx(:,k),oneyx(:,k),'.','color',mymapboy(11,:),'markersize',20)
    plot(twoxx(:,k),twoyx(:,k),'.','color',mymapboy(12,:),'markersize',20)
    plot(threexx(:,k),threeyx(:,k),'.','color',mymapboy(13,:),'markersize',20)
    plot(fourxx(:,k),fouryx(:,k),'.','color',mymapboy(14,:),'markersize',20)
    plot(fivexx(:,k),fiveyx(:,k),'.','color',mymapboy(15,:),'markersize',20)
    plot(sixxx(:,k),sixyx(:,k),'.','color',mymapboy(16,:),'markersize',20)
    plot(sevenxx(:,k),sevenyx(:,k),'.','color',mymapboy(17,:),'markersize',20)
    plot(eightxx(:,k),eightyx(:,k),'.','color',mymapboy(18,:),'markersize',20)
    plot(ninexx(:,k),nineyx(:,k),'.','color',mymapboy(19,:),'markersize',20)
    plot(tenxx(:,k),tenyx(:,k),'.','color',mymapboy(20,:),'markersize',20)
    
    set(gca,'NextPlot','replacechildren')
    frame = getframe(gcf);
    writeVideo(v,frame);
   
    
end

close(v)

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
