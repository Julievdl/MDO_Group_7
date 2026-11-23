clc
clear

%%
Sref = 283.3;
s = 47.57/2;
yk = 7.92;
sweepLE = 35.2;
sweepTE_in = 10.5;
lambda_total = 0.207;

d = yk*(tand(sweepLE)-tand(sweepTE_in));

cr = (Sref + d*s) / ...
     ((1-lambda_total)*yk + (1+lambda_total)*s);

ck = cr-d;
ct = lambda_total*cr;
lambda_out = ct/ck;

S = (cr+ck)*yk + (ck+ct)*(s-yk);

x = [0, yk*tand(sweepLE), s*tand(sweepLE), (s*tand(sweepLE))+ct, (yk*tand(sweepLE))+ck, cr, 0];
y = [0, yk,  s,   s,   yk,  0, 0];

plot(y,x,'o-')
axis equal
set(gca,'YDir','reverse')
xlabel('Spanwise coordinate y [m]')
ylabel('Rearward coordinate x [m]')
grid on

%% MAC calculation

MAC = 2/(3*S) * ( ...
    yk*(cr^2 + cr*ck + ck^2) + ...
    (s-yk)*(ck^2 + ck*ct + ct^2));

MAC_difference_percent = 100*(MAC-6.98)/6.98;

%% Tank Volume and Parameterization
load('whitcomb_CST_baseline.mat','I_tank');
y_end = 0.85*s;
c_end = ck + (ct-ck)*(y_end-yk)/(s-yk);

V_in = I_tank*yk/3*(cr^2 + cr*ck + ck^2);

V_out = I_tank*(y_end-yk)/3 * (ck^2 + ck*c_end + c_end^2);

V_tank = 2*(V_in + V_out);
V_usable = 0.93*V_tank;