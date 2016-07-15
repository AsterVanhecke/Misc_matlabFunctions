

%lg = 1800; % [nm]
dMax = 300; % [nm]
wMax = 0.9; % [adimensional]
Kl = 0.008; % [min]


%-------------------------------------------------------------------------%
% WT

%Kl = 0.014; % [min]
tc = 67; % [min]
tg = 400; % [min]
lc = 2167; % [nm]

t = (tc:1:tg);


l = @(T) exp(Kl*T) + lc + 3*T;

lF = @(T) dMax-dMax*( (tg-T)/(tg -tc));
%l = @(T) 3*T;
%l = @(T) T;
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);
c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));

figure, plot(t, c(t)),
hold on
 plot(t, l(t))
legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t, l(t));
[dldtF, d2ldt2] = dfdxc(t, lF(t));
[dcdtF, d2cdt2] = dfdxc(t, cF(t));
[dcdt, d2cdt2] = dfdxc(t, c(t));

x1 = 1;
y1 = 2;
xF1 = 1;
yF1 = 2;
for i = 2: length(t)
%     xEl= x(i-2)- x(i-1);
%     yEl= y(i-2)- y(i-1);

        xF1(i) = xF1(i-1) + dldtF(i-1) + dldt(i-1);
        yF1(i) = yF1(i-1) + dcdtF(i-1);
        x1(i) = x1(i-1) + dldt(i-1)*(t(i)-t(i-1));
        y1(i) = y1(i-1) + dcdt(i-1);
        
%     if i<=3
%         xF1(i) = xF1(i-1) + dldt(i-1);
%         yF1(i) = yF1(i-1) + dcdtF(i-1);
%         x1(i) = x1(i-1) + dldt(i-1);
%         y1(i) = y1(i-1) + dcdt(i-1);
%     
%     else
% %         xF1(i) = xF1(i-1) + dldt(i-1) + ( xF1(i-2) - xF1(i-1) );
% %         yF1(i) = yF1(i-1) + dcdtF(i-1) + ( yF1(i-2) - yF1(i-1) );
% %         x1(i) = x1(i-1) + dldt(i-1) + ( x1(i-2) - x1(i-1) );
% %         y1(i) = y1(i-1) + dcdt(i-1) + ( y1(i-2) - y1(i-1) );
%         
%         xF1(i) = xF1(i-1) + dldt(i-1) ;
%         yF1(i) = yF1(i-1) + dcdtF(i-1)+ dldt(i-1) ;
%         x1(i) = x1(i-1) + dldt(i-1) + dldt(i-1);
%         y1(i) = y1(i-1) + dcdt(i-1)+ dldt(i-1);
%     end
end

figure, plot(x1, y1, '*');
legend('bacteria contour')
% figure, plot(t, dcdt,'*')
legend ('constrictio derivative')

%-----------------------------------------------------------------------%
% 3 MUT

%Kl = 0.014;
tc = 95;
tg = 253;
lc = 2395; % [nm]

t = (tc:1:tg);

l = @(T) exp(Kl*T) + lc +3*T;
lF = @(T) dMax- dMax*( (tg-T)/(tg -tc));

%l = @(T) 3*T;
%l = @(T) T;

c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);

figure, plot(t, c(t)),
hold on
plot(t, l(t))
legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t, l(t));
[dldtF, d2ldt2] = dfdxc(t, lF(t));
[dcdt, d2cdt2] = dfdxc(t, c(t));
[dcdtF, d2cdt2] = dfdxc(t, cF(t));

x4=1;
y4=2;
xF4=1;
yF4=2;

for i =2: length(t)
%     xEl= x(i-2)- x(i-1);
%     yEl= y(i-2)- y(i-1);
    x4(i) = x4(i-1) + dldt(i-1)*(t(i)-t(i-1));
    y4(i) = y4(i-1) + dcdt(i-1);
    xF4(i) = xF4(i-1) + dldtF(i-1)+ dldt(i-1);
    yF4(i) = yF4(i-1) + dcdtF(i-1);
end

% figure, plot(x4, y4, '*');
% legend('bacteria contour')
% figure, plot(t, dcdt,'*')
% legend ('constrictio derivative')

%-----------------------------------------------------------------------%
% MUT FtsW

%Kl = 0.014;
tc = 148;
tg = 360;
lc = 2180; % [nm]

t1 = (tc:1:tg);


l = @(T) exp(Kl*T) + lc + 3*T;
lF = @(T)  dMax-dMax*( (tg-T)/(tg -tc));

%l = @(T) 3*T;
%l = @(T) T;
c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2 );
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));

figure, plot(t1, c(t1)),
hold on
 plot(t1, l(t1))
legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t1, l(t1));
[dldtF, d2ldt2] = dfdxc(t1, lF(t1));
[dcdt, d2cdt2] = dfdxc(t1, c(t1));
[dcdtF, d2cdt2] = dfdxc(t1, cF(t1));

x2=1;
y2=2;
xF2=1;
yF2=2;
for i = 2: length(t1)
%     xEl= x(i-2)- x(i-1);
%     yEl= y(i-2)- y(i-1);
    x2(i) = x2(i-1) + dldt(i-1)*(t1(i)-t1(i-1));
    y2(i) = y2(i-1) + dcdt(i-1);
    xF2(i) = xF2(i-1) + dldtF(i-1)+ dldt(i-1);
    yF2(i) = yF2(i-1) + dcdtF(i-1);
end

% figure, plot(t1, dcdt,'*')
% legend ('constrictio derivative')
% 
% figure, plot(x1, y1, 'b*');
% hold on
%  plot(x2, y2, 'r*');
% legend('bacteria contour tc=67 and tg=400', 'bacteria contour tc=93 and tg=290')


%-----------------------------------------------------------------------%
% FtsI

% Kl = 0.014;
tc = 137;
tg = 337;
lc = 2343; % [nm]

t2 = (tc:1:tg);


l = @(T) exp(Kl*T) + lc +3*T;
lF = @(T)  dMax-dMax*( (tg-T)/(tg -tc));

%%l = @(T) T;
c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)) );
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2 );

figure, plot(t2, c(t2)),
hold on
 plot(t2, l(t2))
legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t2, l(t2));
[dldtF, d2ldt2] = dfdxc(t2, lF(t2));
[dcdt, d2cdt2] = dfdxc(t2, c(t2));
[dcdtF, d2cdt2] = dfdxc(t2, cF(t2));

x3=1;
y3=2;
xF3=1;
yF3=2;
% x(2)=2;
% y(2)=1;
for i =2: length(t2)
%     xEl= x(i-2)- x(i-1);
%     yEl= y(i-2)- y(i-1);
    x3(i) = x3(i-1) + dldt(i-1)*(t2(i)-t2(i-1));
    y3(i) = y3(i-1) + dcdt(i-1);
    xF3(i) = xF3(i-1) + dldtF(i-1)+ dldt(i-1);
    yF3(i) = yF3(i-1) + dcdtF(i-1);
end

figure, plot(t2, dcdt,'*')
legend ('constrictio derivative')

figure, plot(x1, y1, 'b*'),
hold on
plot(x1, 2*y1(end)- y1, 'b*'),
hold all
plot(x2, y2, 'r*');
plot(x2, 2*y2(end)- y2  , 'r*'),
plot(x3, y3, 'g*');
plot(x3, 2*y3(end)- y3, 'g*'),
plot(x4, y4, 'k*');
plot(x4, 2*y4(end)- y4, 'k*'),
xlabel('[nm]')
ylabel('[nm]')
axis equal
title('2D area model')
legend('bacteria contour tc=67 and tg=400, WT', 'bacteria contour tc=67 and tg=400, WT',...
    'bacteria contour tc=148 and tg=360, FtsW(A246T)', 'bacteria contour tc=148 and tg=360, FtsW(A246T)',...
'bacteria contour tc=137 and tg=337 (FtsI I45V)', 'bacteria contour tc=137 and tg=337 (FtsI I45V)',...
'bacteria contour tc=95 and tg=253 FtsW39', 'bacteria contour tc=95 and tg=253 FtsW39')

figure, plot(xF1, yF1, 'b*'),
hold on
plot(xF1, 2*yF1(end)- yF1, 'b*'),
hold all
plot(xF2, yF2, 'r*');
plot(xF2, 2*yF2(end)- yF2  , 'r*'),
plot(xF3, yF3, 'g*');
plot(xF3, 2*yF3(end)- yF3, 'g*'),
plot(xF4, yF4, 'k*');
plot(xF4, 2*yF4(end)- yF4, 'k*'),
xlabel('[nm]')
ylabel('[nm]')
axis equal
title('Feingold Model')
legend('bacteria contour tc=67 and tg=400, WT', 'bacteria contour tc=67 and tg=400, WT',...
    'bacteria contour tc=148 and tg=360, FtsW(A246T)', 'bacteria contour tc=148 and tg=360, FtsW(A246T)',...
'bacteria contour tc=137 and tg=337 (FtsI I45V)', 'bacteria contour tc=137 and tg=337 (FtsI I45V)',...
'bacteria contour tc=95 and tg=253 FtsW39', 'bacteria contour tc=95 and tg=253 FtsW39')

