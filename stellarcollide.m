function stellarcollide()
%Begin by randomizing mass and location of planets for star 1
nPlanets1 = 8; nPlanets2 = 8;
nTotBodies = nPlanets1+nPlanets2+2;
G = 6.67384e-20; %Gravitational constant in km^3 kg^-1 s^-2
nex1 = 2e7*2;
nex2 = 2e7*2;
MaxDistRatio=4;
Thermal_energy_system = 0;
k = 0;
mmin = 1.3e22;
mmax = 1.9e27*10;
m1 = random('uniform',mmin,mmax,[1,nPlanets1]);
x1 = random('uniform',-nex1,nex1,[1,nPlanets1]);
y1 = random('uniform',-nex1,nex1,[1,nPlanets1]);
nIter = 2000;
dt = 2000;

ms1 = sum(m1)/.01;
%same for star 2
m2 = random('uniform',mmin,mmax,[1,nPlanets2]);
x2 = random('uniform',-nex2,nex2,[1,nPlanets2]);
y2 = random('uniform',-nex2,nex2,[1,nPlanets2]);

ms2 = sum(m2)/.01;
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
rhat = zeros(2,nTotBodies);
v = zeros(1,nTotBodies);
vx = zeros(1,nTotBodies+1);
vy = zeros(1,nTotBodies+1);
%Calculate velocity according to F_c = mv^2/r
for i = 1:nTotBodies
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
a = [0 -1; 1 0]; %Matrix to make perpendicular vector
vhat = a*rhat;
for i = 1:2:nTotBodies
    vhat(1,i) = -vhat(1,i);
    vhat(2,i) = -vhat(2,i);
end
for i = 1:nTotBodies
    vx(i) = v(i)*vhat(1,i);
    vy(i) = v(i)*vhat(2,i);
end
R = max(rs)*MaxDistRatio;
cm = (ms2 +sum(m2))*R/(ms1+ms2+sum(m1)+sum(m2));
for i = 1:nTotBodies
    if i > nPlanets1+1
        x_comb(i) = x_comb(i) + (R-cm);
    else
        x_comb(i) = x_comb(i) - cm;
    end
end
a_Gx = zeros(1,nTotBodies);
a_Gy = zeros(1,nTotBodies);
for i = 1:nTotBodies
    for j = i+1:nTotBodies
        r(i,j) = sqrt((x_comb(j)-x_comb(i))^2 + (y_comb(j)-y_comb(i))^2);
        r(j,i) = r(i,j);
        rx(i,j) = (x_comb(i) - x_comb(j))/r(i,j);
        rx(j,i) = -rx(i,j);
        ry(i,j) = (y_comb(i) - y_comb(j))/r(i,j);
        ry(j,i) = -ry(i,j);
    end
    for j = 1:nTotBodies
        if i ~= j
            a_Gx(i) = a_Gx(i) - G*(m_comb(j)/(r(i,j)^2))*rx(i,j);
            a_Gy(i) = a_Gy(i) - G*(m_comb(j)/(r(i,j)^2))*ry(i,j);
        end
    end
end

for s = 1:nIter
    for i = 1:nTotBodies
        x_comb(i) = x_comb(i) + dt*vx(i) + (dt^2)*a_Gx(i)/2;
        y_comb(i) = y_comb(i) + dt*vy(i) + (dt^2)*a_Gy(i)/2;
    end
    for i = 1:nTotBodies
        for j = i+1:nTotBodies
            r(i,j) = sqrt((x_comb(j)-x_comb(i))^2 + (y_comb(j)-y_comb(i))^2);
            r(j,i) = r(i,j);
            rx(i,j) = (x_comb(i) - x_comb(j))/r(i,j);
            rx(j,i) = -rx(i,j);
            ry(i,j) = (y_comb(i) - y_comb(j))/r(i,j);
            ry(j,i) = -ry(i,j);
        end
    end
    
    if r(nTotBodies,nPlanets1+1) < 5000000
        fprintf('Yes')
        m_new = m_comb(nTotBodies) + m_comb(nPlanets1+1);
        Thermal_energy_system = .5*m_comb(nTotBodies)*(vx(nTotBodies)^2 + vy(nTotBodies)^2)/2 + .5*m_comb(nPlanets1+1)*(vx(nPlanets1+1)^2 + vy(nPlanets1+1)^2)/2;
        %         vx(nTotBodies+1) = m_comb(nTotBodies)*vx(nTotBodies)/m_new + m_comb(nPlanets1+1)*vx(nPlanets1+1)/m_new;
        %         vy(nTotBodies+1) = m_comb(nTotBodies)*vy(nTotBodies)/m_new + m_comb(nPlanets1+1)*vy(nPlanets1+1)/m_new;
        vx(nTotBodies) = 0;
        vy(nTotBodies) = 0;
        k = 2;
        m_comb(nTotBodies) = m_new; m_comb(nPlanets1+1) = 0; vx(nPlanets1+1) = 0; vy(nPlanets1+1) = 0; a_Gx(nPlanets1+1) = 0; a_Gy(nPlanets1+1) = 0; x_comb(nPlanets1+1) = 6e15;
    end
    for j = 1:nTotBodies
        a_Gx_t(j) = a_Gx(j);
        a_Gy_t(j) = a_Gy(j);
    end
    a_Gx = zeros(1,nTotBodies);
    a_Gy = zeros(1,nTotBodies);
    for i = 1:nTotBodies
        for j = 1:nTotBodies
            if i ~= j
                a_Gx(i) = a_Gx(i) - G*(m_comb(j)/(r(i,j)^2))*rx(i,j);
                a_Gy(i) = a_Gy(i) - G*(m_comb(j)/(r(i,j)^2))*ry(i,j);
            end
        end
    end
    for i = 1:nTotBodies
        if k ~= 0
            vx(i) = vx(i) + dt*a_Gx(i);
            vy(i) = vy(i) + dt*a_Gy(i);
            vx(nTotBodies) = 0;
            vy(nTotBodies) = 0;
            k = k-1;
        else
            vx(i) = vx(i) + dt*(a_Gx(i)+a_Gx_t(i))/2;
            vy(i) = vy(i) + dt*(a_Gy(i)+a_Gy_t(i))/2;
            if vx(nTotBodies+1) ~= 0
                vx(nTotBodies) = vx(nTotBodies+1);
                vx(nTotBodies+1) = 1;
                vy(nTotBodies) = vy(nTotBodies+1);
                vy(nTotBodies+1) = 0;
            end
        end
    end
    
    clf
    plot(x_comb(1:nPlanets1),y_comb(1:nPlanets1),'+')
    hold on
    plot(x_comb(nPlanets1+1),y_comb(nPlanets1+1),'r+')
    plot(x_comb(nPlanets1+2:nTotBodies-1),y_comb(nPlanets1+2:nTotBodies-1),'o')
    plot(x_comb(nTotBodies),y_comb(nTotBodies),'or')
    axis([-cm-max(rs)*2 (R-cm)+max(rs)*2 -max(rs)*1.5*3 max(rs)*1.5*2])
    movie(s) = getframe;
    
end