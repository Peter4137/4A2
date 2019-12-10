close all;
clear
fin1 = fopen('var_dt/euler.mat','r');
data1 = fscanf(fin1,'%f');
ni1 = data1(1,1);
nj1 = data1(2,1);
idone = 2;
x1 = zeros(ni,nj);
y1 = zeros(ni,nj);
ro1 = zeros(ni,nj);
vx1 = zeros(ni,nj);
vy1 = zeros(ni,nj);
p1 = zeros(ni,nj);
%
for i=1:ni
  for j=1:nj
    x1(i,j)  = data1(idone+1,1);
    y1(i,j)  = data1(idone+2,1);
    ro1(i,j) = data1(idone+3,1);
    vx1(i,j) = data1(idone+4,1);
    vy1(i,j) = data1(idone+5,1);
    p1(i,j)  = data1(idone+6,1);
    idone = idone+6;
  end 
end
fclose(fin1);

fin2 = fopen('var_dt/euler.mat','r');
data2 = fscanf(fin,'%f');
ni2 = data2(1,1);
nj2 = data2(2,1);
idone = 2
x2 = zeros(ni,nj);
y2 = zeros(ni,nj);
ro2 = zeros(ni,nj);
vx2 = zeros(ni,nj);
vy2 = zeros(ni,nj);
p2 = zeros(ni,nj);
%
for i=1:ni
  for j=1:nj
    x2(i,j)  = data2(idone+1,1);
    y2(i,j)  = data2(idone+2,1);
    ro2(i,j) = data2(idone+3,1);
    vx2(i,j) = data2(idone+4,1);
    vy2(i,j) = data2(idone+5,1);
    p2(i,j)  = data2(idone+6,1);
    idone = idone+6;
  end 
end
fclose(fin2);

plot(vx1(ni1,:)-vx2(ni2,:), y(ni1,:))
pause



% Contour Static Pressure
%
figure('Name','Static Pressure');
hold on
daspect([1 1 1]);
plot(x(:,1),y(:,1),'b','Linewidth',1);
plot(x(:,nj),y(:,nj),'b','LineWidth',1);
pmax = round(max(max(p)),1,'significant');
pmin = round(min(min(p)),1,'significant');
delp = (pmax-pmin)/pmin;
V = linspace(pmin,pmax,21);
[C,h] = contourf(x,y,p,V,'LineWidth',2);
%clabel(C,'manual','FontSize',14);
hold off;
pause
%
% Plot grid
%
figure('Name','Grid');
hold on;
daspect([1 1 1]);
for j=1:nj
  plot(x(:,j),y(:,j),'b');
end
for i=1:ni
  plot(x(i,:),y(i,:),'b');
end
hold off;
%
% Contour derived variables
%
gam = 1.4;
a = sqrt(gam*p(1,1)/ro(1,1));
p0 = zeros(ni,nj);
mach = zeros(ni,nj);
for i=1:ni
  for j=1:nj
    mach(i,j) = sqrt(ro(i,j)*(vx(i,j)^2+vy(i,j)^2)/(gam*p(i,j)));
    p0(i,j) = p(i,j)*(1.+.5*(gam-1.)*mach(i,j)^2)^(gam/(gam-1));
  end 
end
%
figure('Name','Stagnation Pressure');
hold on
daspect([1 1 1]);
contourf(x,y,p0,20,'LineWidth',2);
plot(x(:,1),y(:,1),'b','Linewidth',1);
plot(x(:,nj),y(:,nj),'b','LineWidth',1);
p0max = max(max(p0))
p0min = min(min(p0))
%
figure('Name','Mach Number');
hold on
daspect([1 1 1]);
contourf(x,y,mach,20,'LineWidth',2);
plot(x(:,1),y(:,1),'b','Linewidth',1);
plot(x(:,nj),y(:,nj),'b','LineWidth',1);
machmax = max(max(mach))
machmin = min(min(mach))