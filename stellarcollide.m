%This is a graphical simulation of the collision of two stellar systems

function stellarcollide()
nPlanets1 = 50; nPlanets2 = 80; %The number of planets in each stellar system
nTotBodies = nPlanets1+nPlanets2+2; %total number of gravitational bodies
G = 6.67384e-20; %Gravitational constant in km^3 kg^-1 s^-2
MaxDistRatio=4; %A variable for controlling the starting distance between the two suns
StarMassRatio=5000; %the ratio between a stars mass and the sum of its planet's masses
PlanetRadiusRatio = 2; %ratio of stellar system radius ranges
nex1 = 4e7; %these two variables control the range of radius allowed between the stars and their planets
nex2 = nex1*PlanetRadiusRatio;
mmin = 1.3e22; %these two variables control the range of masses for the planets
mmax = 1.9e27*10;
m1 = random('uniform',mmin,mmax,[1,nPlanets1]); %randomly assign masses and radii for first stellar system
x1 = random('uniform',-nex1,nex1,[1,nPlanets1]);
y1 = random('uniform',-nex1,nex1,[1,nPlanets1]);
nIter = 500; %number of iterations to perform
dt = 2000; %arbitary

ms1 = sum(m1)*StarMassRatio;
%same for star 2
m2 = random('uniform',mmin,mmax,[1,nPlanets2]); %randomly assign planet masses and radii for first stellar system
x2 = random('uniform',-nex2,nex2,[1,nPlanets2]);
y2 = random('uniform',-nex2,nex2,[1,nPlanets2]);

ms2 = sum(m2)*StarMassRatio;

%%%%%% this code takes the now generated stellar system data and combines them into arrays containing     
     % all gravitational bodies in order to easily loop over them
x_comb(1:nPlanets1) = x1;
x_comb(nPlanets1+1) = 0;
x_comb(nTotBodies) = 0;
x_comb(nPlanets1+2:(nTotBodies-1)) = x2;
y_comb(1:nPlanets1) = y1;
y_comb(nPlanets1+1) = 0;
y_comb(nPlanets1+2:(nTotBodies-1)) = y2;
y_comb(nTotBodies) = 0;
m_comb(1:nPlanets1) = m1;
m_comb(nPlanets1+1) = ms1;
m_comb(nPlanets1+2:nTotBodies-1) = m2;
m_comb(nTotBodies) = ms2;
%%%%%
rhat = zeros(2,nTotBodies); %contains initial velocity unit vectors
v = zeros(1,nTotBodies);
vx = zeros(1,nTotBodies+1);
vy = zeros(1,nTotBodies+1);

for i = 1:nTotBodies %This for loop calculates velocities of the planets according to F_c = mv^2/r 
                     %so as to be in stable orbit when simulation begins
    if i ~= nPlanets1+1
        if i ~= nTotBodies
            rs(i) = sqrt(y_comb(i)^2 + x_comb(i)^2);
            if (i > nPlanets1)
                v(i) = sqrt(G*ms2/rs(i));
            else
                v(i) = sqrt(G*ms1/rs(i)); 
            end
            rhat(1,i) = -x_comb(i)/rs(i);
            rhat(2,i) = -y_comb(i)/rs(i);
        end
    end
end
a = [0 -1; 1 0];
vhat = a*rhat; %creates velocity unit vectors perpendicular to radius
for i = 1:2:nTotBodies %make every other orbit in opposite direction (just for fun, this almost never happens in reality)
    vhat(1,i) = -vhat(1,i);
    vhat(2,i) = -vhat(2,i);
end
for i = 1:nTotBodies
    vx(i) = v(i)*vhat(1,i);
    vy(i) = v(i)*vhat(2,i);
end
R = max(rs)*MaxDistRatio; %R holds the distance between the two stars
cm = (ms2 +sum(m2))*R/(ms1+ms2+sum(m1)+sum(m2)); %calculate center of mass of system
for i = 1:nTotBodies
    if i > nPlanets1+1
        x_comb(i) = x_comb(i) + (R-cm); %shift planets in the x over into their places around the stars. 
                                        %until now all position data for planets was reletive to their star, not absolute.
    else
        x_comb(i) = x_comb(i) - cm;
    end
end
a_Gx = zeros(1,nTotBodies); %accelaration due to gravity
a_Gy = zeros(1,nTotBodies);
for i = 1:nTotBodies
    for j = i+1:nTotBodies %this loop calculates the initial reletive distances between all planets
        r(i,j) = sqrt((x_comb(j)-x_comb(i))^2 + (y_comb(j)-y_comb(i))^2);
        r(j,i) = r(i,j);
        rx(i,j) = (x_comb(i) - x_comb(j))/r(i,j);
        rx(j,i) = -rx(i,j);
        ry(i,j) = (y_comb(i) - y_comb(j))/r(i,j);
        ry(j,i) = -ry(i,j);
    end
    for j = 1:nTotBodies %this loop calculates the initial accelaration due to gravity for all planets
        if i ~= j
            a_Gx(i) = G*(m_comb(j)/(r(i,j)^2))*rx(i,j);
            a_Gy(i) = G*(m_comb(j)/(r(i,j)^2))*ry(i,j);
        end
    end
end

for s = 1:nIter %this is the main loop and propagates the system in time after the previous set-up phase
    for i = 1:nTotBodies %the Verlet algorithm for propagation of the gravitational bodies
        x_comb(i) = x_comb(i) + dt*vx(i) + (dt^2)*a_Gx(i)/2;
        y_comb(i) = y_comb(i) + dt*vy(i) + (dt^2)*a_Gy(i)/2;
    end
    for i = 1:nTotBodies %recalculate new relative distances between gravitational bodies
        for j = i+1:nTotBodies
            r(i,j) = sqrt((x_comb(j)-x_comb(i))^2 + (y_comb(j)-y_comb(i))^2);
            r(j,i) = r(i,j);
            rx(i,j) = (x_comb(i) - x_comb(j))/r(i,j);
            rx(j,i) = -rx(i,j);
            ry(i,j) = (y_comb(i) - y_comb(j))/r(i,j);
            ry(j,i) = -ry(i,j);
        end
    end
    for j = 1:nTotBodies
        a_Gx_t(j) = a_Gx(j); %save accelaration from previous iteration 
        a_Gy_t(j) = a_Gy(j);
    end
    a_Gx = zeros(1,nTotBodies); %get rid of accelaration data from the previous timestep
    a_Gy = zeros(1,nTotBodies);
    for i = 1:nTotBodies
        for j = 1:nTotBodies
            if i ~= j
                a_Gx(i) = a_Gx(i) - G*(m_comb(j)/(r(i,j)^2))*rx(i,j); %%calculate new accelaration according to ay_k+1 = ay_k - GmM/((r^2)*ry and 
                a_Gy(i) = a_Gy(i) - G*(m_comb(j)/(r(i,j)^2))*ry(i,j); %corresponding equations for x
            end
        end
    end
    for i = 1:nTotBodies
            vx(i) = vx(i) + dt*(a_Gx(i)+a_Gx_t(i))/2; %verlet algorithm for velocity, using the average between k and k+1 accelaration
            vy(i) = vy(i) + dt*(a_Gy(i)+a_Gy_t(i))/2;
            if vx(nTotBodies+1) ~= 0
                vx(nTotBodies) = vx(nTotBodies+1);
                vx(nTotBodies+1) = 1;
                vy(nTotBodies) = vy(nTotBodies+1);
                vy(nTotBodies+1) = 0;
            end
    end
    
    clf
    plot(x_comb(1:nPlanets1),y_comb(1:nPlanets1),'+')
    hold on
    plot(x_comb(nPlanets1+1),y_comb(nPlanets1+1),'r+')
    plot(x_comb(nPlanets1+2:nTotBodies-1),y_comb(nPlanets1+2:nTotBodies-1),'o')
    plot(x_comb(nTotBodies),y_comb(nTotBodies),'or')
    axis([-cm-max(rs)*2 (R-cm)+max(rs)*2 -max(rs)*4.5 max(rs)*3]) %adjusts the axis of the plot dynamically so everything is always visible
    movie(s) = getframe; %This causes the plot to be displayed every iteration instead of just at the end
    
end
