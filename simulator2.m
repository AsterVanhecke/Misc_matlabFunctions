% ------------------------------------------------------------------------%
%                 POLE SHAPE SIMULATOR
%-------------------------------------------------------------------------%
% This function simulates how the new two poles of a Caulobacter 
% are built at the bacteria septum when the constriction starts.
% So, the starting bacteria length is the length at Tc, at the time when
% the constriction starts.
% The main hypotesis is that when the constriction starts all the
% constribution to the elongation is given at the septum.
%
% Here, two models are taken into account: 
% - a 2D model where the first hypotesis is that the elogation contribution 
%   during the constriction is entirely given by the septum site. 
%   The elongation is an exponential process that is due to an insertion of 
%   PG everywhere along the pole that is builded up during the contriction 
% - Feingold model where the "pole" elongation is only due to the fact that 
%   the septum surface is growing in 3D. (See definition of h(t))
%
% This simulator builds up the pole contour translating an initial point
% using a costriction vector and an elongation vector. 
% The model takes into account that, in the constriction site, the PG 
% insertion occurs at each time step not only where a new point is created 
% but also where they were already created during the constriction.
%
% Instead of to create at each step a new point and, simultaneuly, to add 
% an elongation vector to all the ones already created we are creating the 
% new point includeing immeditatly all the elongation contributions it had 
% during each step of the contriction.
%
% - the constriction vector has a modulus depending on the rate of the
%   constriction function and a direction perpendicular to the bacteria lenght
% - the elongation vector has a modulus depending on the total number of steps 
%   have to be done to finish the constriction. The length of these steps
%   depends on the elongation step. The direction of this vector is
%   parallel to the lenght of the bacteria
%
% The cell elongation is considered exponential.

 

%
%-------------------------------------------------------------------------%

% General pramters definition

%lc = 1800; % [nm]
dMax = 240; % [nm] divided by 2
wMax = 0.9; % [adimensional]
%Kl = 0.008; % [min]
Kl = 0.01; % [min]

%-------------------------------------------------------------------------%
% WT
%-------------------------------------------------------------------------%

% General WT parameters
%Kl = 0.014; % [min]
tc = 67; % [min]
tg = 400; % [min]
lc = 2167; % [nm]

% Define the time vector
t = (tc:1:tg);

% General exponential elongation 
% We are looking to only one pole --> elongation / 2
l = @(T) ((exp(Kl*T)) )/2;

% General linear elongation 
%l = @(T) Kl*T + lc;

% Septum Contribution to the Elongation (Feingold model - h(t))
lF = @(T) (dMax-dMax*( (tg-T)/(tg -tc)));

% Septum constriction (Feingold model and 2D area model)
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);
c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));

% DEBUG
% figure, plot(t, c(t)),
% hold on
% plot(t, l(t))
% legend('constriction vs time', 'length vs time')

% Derivative(rate) elongation and constriction (Feingold model and 2D area model)
[dldt, d2ldt2] = dfdxc(t, l(t));
[dldtF, d2ldt2] = dfdxc(t, lF(t));
[dcdtF, d2cdt2] = dfdxc(t, cF(t));
[dcdt, d2cdt2] = dfdxc(t, c(t));

% Choose an initial point
x1 = 1;
y1 = 2;
xF1 = 1;
yF1 = 2;

% Build up the pole contour adding a costriction vector and an elongation
% vector depending on the PG insertion along the bacteria length

% Initialize lTot
lTot = 0;
lTotF = 0;
for i = 2: length(t) % cycle over the time steps

    % Cycle to grab the right elongation at each step 
    % that has to be summed up to this new point
    for dlIdx = i: length(t) % from the elongation step i to the last
        % 2D model
        lTot = dldt(i-1)./i + lTot;
        
        % Feingol model
        %lTotF = dldtF(i-1)./i + lTotF;
        %lTotF = (dldtF(i-1) + dldt(i-1))./i + lTotF;
    end
    
    % Feingold model
    xF1(i) = xF1(i-1) + dldtF(i-1)*(t(i)-t(i-1));
    yF1(i) = yF1(i-1) + dcdtF(i-1);
    % 2D model
    x1(i) = x1(i-1) + lTot*(t(i)-t(i-1));
    y1(i) = y1(i-1) + dcdt(i-1);
end % end cycle over the time step

% contour top (left and right)
x1TempTop = [x1'; 2*x1(end) - flipud(x1')];
y1TempTop = [y1'; flipud(y1')];
% contour bottom (left and right)
x1TempBot = [x1'; 2*x1(end) - flipud(x1')];
y1TempBot = [2*y1(end) - y1'; flipud(2*y1(end) - y1')];

% contour top (left and right)
x1TempTopF = [xF1'; 2*xF1(end) - flipud(xF1')];
y1TempTopF = [yF1'; flipud(yF1')];
% contour bottom (left and right)
x1TempBotF = [xF1'; 2*xF1(end) - flipud(xF1')];
y1TempBotF = [2*yF1(end) - yF1'; flipud(2*yF1(end) - yF1')];

% contour with both right and left pole 
contourWTAll = [x1TempTop, y1TempTop, x1TempBot, y1TempBot];
contourWTFAll = [x1TempTopF, y1TempTopF, x1TempBotF, y1TempBotF];

% contour of only the left pole
contourWTHalf = [x1', y1', x1', 2*y1(end) - y1'];
contourWTFHalf = [xF1', yF1', xF1', 2*yF1(end) - yF1'];

%-------------------------------------------------------------------------%
% Dynamic pole WT
%-------------------------------------------------------------------------%
% figure,
% col=hsv(length(t));
% for timeIdx = 3:1:length(t)
%     
%     % Define the time vector
%         tTemp = t(1:timeIdx);
% 
%         % General exponential elongation 
%         % We are looking to only one pole --> elongation / 2
%         l = @(T) ((exp(Kl*T)) + lc)/2;
% 
%         % General linear elongation 
%         %l = @(T) Kl*T + lc;
% 
%         % Septum Contribution to the Elongation (Feingold model - h(t))
%         lF = @(T) (dMax-dMax*( (tg-T)/(tg -tc)));
% 
%         % Septum constriction (Feingold model and 2D area model)
%         cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);
%         c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));
% 
%         % Derivative(rate) elongation and constriction (Feingold model and 2D area model)
%         [dldt, d2ldt2] = dfdxc(tTemp, l(tTemp));
%         [dldtF, d2ldt2] = dfdxc(tTemp, lF(tTemp));
%         [dcdtF, d2cdt2] = dfdxc(tTemp, cF(tTemp));
%         [dcdt, d2cdt2] = dfdxc(tTemp, c(tTemp));
% 
%         % Choose an initial point
%         xd1 = 1;
%         yd1 = 2;
%         xdF1 = 1;
%         ydF1 = 2;
% 
%         % Build up the pole contour adding a costriction vector and an elongation
%         % vector depending on the PG insertion along the bacteria length
% 
%         % Initialize lTot
%         lTot = 0;
%         lTotF = 0;
%         for i = 2: length(tTemp) % cycle over the time steps
% 
%             % Cycle to grab the right elongation at each step 
%             % that has to be summed up to this new point
%             for dlIdx = i: length(tTemp) % from the elongation step i to the last
%                 % 2D model
%                 lTot = dldt(i-1)./i + lTot;
% 
%                 % Feingol model
%                 %lTotF = dldtF(i-1)./i + lTotF;
%                 %lTotF = (dldtF(i-1) + dldt(i-1))./i + lTotF;
%             end
% 
%             % Feingold model
%             xdF1(i) = xdF1(i-1) + dldtF(i-1)*(tTemp(i)-tTemp(i-1));
%             ydF1(i) = ydF1(i-1) + dcdtF(i-1);
%             % 2D model
%             xd1(i) = xd1(i-1) + lTot*(tTemp(i)-tTemp(i-1));
%             yd1(i) = yd1(i-1) + dcdt(i-1);
%         end % end cycle over the time step
% 
%         % contour top (left and right)
%         xd1TempTop = [xd1' - xd1(end) + x1(end); 2*xd1(end) - flipud(xd1') - xd1(end) + x1(end)];
%         yd1TempTop = [yd1'; flipud(yd1')];
%         % contour bottom (left and right)
%         xd1TempBot = [xd1' - xd1(end) + x1(end); 2*xd1(end) - flipud(xd1') - xd1(end) + x1(end)];
%         yd1TempBot = [-2*dMax - yd1'; flipud(-2*dMax - yd1')];
% 
%         % contour top (left and right)
%         xd1TempTopF = [xdF1' - xdF1(end) + xF1(end); 2*xdF1(end) - flipud(xdF1') + xF1(end)];
%         yd1TempTopF = [ydF1'; flipud(ydF1')];
%         % contour bottom (left and right)
%         xd1TempBotF = [xdF1' - xdF1(end) + xF1(end); 2*xdF1(end) - flipud(xdF1') + xF1(end)];
%         yd1TempBotF = [-2*dMax - ydF1'; flipud(-2*dMax - ydF1')];
% 
%         % contour with both right and left pole 
%         contourWTAlld = [xd1TempTop, yd1TempTop, xd1TempBot, yd1TempBot];
%         contourWTFAlld = [xd1TempTopF, yd1TempTopF, xd1TempBotF, yd1TempBotF];
% 
%         % contour of only the left pole
%         contourWTHalfd = [xd1' - xd1(end) + x1(end), yd1', xd1' - xd1(end) + x1(end),- 2*dMax - yd1'];
%         contourWTFHalfd = [xdF1'- xdF1(end) + xF1(end), ydF1', xdF1' - xdF1(end) + xF1(end),- 2*dMax - ydF1'];
% 
%         % DEBUG
%         subplot(1, 2, 1)
%         plot(contourWTAlld(:, 1), contourWTAlld(:, 2), '*', 'Color',col(i,:), 'markersize', 1),
%         hold on,
%         plot(contourWTAlld(:, 3), contourWTAlld(:, 4), '*', 'Color',col(i,:), 'markersize', 1),
%         axis equal
%         xlabel('length [nm]'),
%         ylabel('length [nm]'),
%         legend('WT simulation'),
%         drawnow
% end % end dynamic pole
%-------------------------------------------------------------------------%

% DEBUG
% figure, plot(x1, y1, '*');
% legend('bacteria contour')
% % figure, plot(t, dcdt,'*')
% legend ('constrictio derivative')

%-------------------------------------------------------------------------%
% General pramters definition

%lc = 1800; % [nm]
dMax = 225; % [nm]
wMax = 0.9; % [adimensional]
%Kl = 0.008; % [min]
Kl = 0.003; % [min]
%-------------------------------------------------------------------------%

lTot = 0;
lTotF = 0;
%-------------------------------------------------------------------------%
% 3 MUT
%-------------------------------------------------------------------------%

%Kl = 0.014;
tc = 118;
tg = 269;
lc = 2395; % [nm]

t = (tc:1:tg);

l = @(T) ((exp(Kl*T)) + lc)/2;
lF = @(T) (dMax-dMax*( (tg-T)/(tg -tc)));

c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);

% DEBUG
% figure, plot(t, c(t)),
% hold on
% plot(t, l(t))
% legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t, l(t));
[dldtF, d2ldt2] = dfdxc(t, lF(t));
[dcdt, d2cdt2] = dfdxc(t, c(t));
[dcdtF, d2cdt2] = dfdxc(t, cF(t));

x4=1;
y4=2;
xF4=1;
yF4=2;

for i =2: length(t)

    for dlIdx = i: length(t) % from the elongation step i to the last
        % 2D model
        lTot = dldt(i-1)./i + lTot;
        
        % Feingol model
        %lTotF = dldtF(i-1)./i + lTotF; 
    end
    x4(i) = x4(i-1) + lTot*(t(i)-t(i-1));
    y4(i) = y4(i-1) + dcdt(i-1);
    
    xF4(i) = xF4(i-1) + dldtF(i-1)*(t(i)-t(i-1));
    yF4(i) = yF4(i-1) + dcdtF(i-1);
end

% contour top (left and right)
x4TempTop = [x4'; x4' + x4(end)];
y4TempTop = [y4'; flipud(y4')];
% contour bottom (left and right)
x4TempBot = [x4'; x4' + x4(end)];
y4TempBot = [2*y4(end) - y4'; flipud(2*y4(end) - y4')];

% contour top (left and right)
x4TempTopF = [xF4'; xF4' + xF4(end)];
y4TempTopF = [yF4'; flipud(yF4')];
% contour bottom (left and right)
x4TempBotF = [xF4'; xF4' + xF4(end)];
y4TempBotF = [2*yF4(end) - yF4'; flipud(2*yF4(end) - yF4')];

% center with the MUT contour (left and right poles)
contourMUTWIAll = [x4TempTop, y4TempTop + (y1(end) - y4(end)), x4TempBot, y4TempBot + (y1(end) - y4(end))];
contourMUTWIFAll = [x4TempTopF, y4TempTopF + (yF1(end) - yF4(end)), x4TempBotF, y4TempBotF + (yF1(end) - yF4(end))];

% contour without rigth part
contourMUTWIHalf = [x4', y4' + (y1(end) - y4(end)), x4', 2*y4(end) - y4' + (y1(end) - y4(end))];
contourMUTWIFHalf = [xF4', yF4' + (yF1(end) - yF4(end)), xF4', 2*yF4(end) - yF4' + (yF1(end) - yF4(end))];

%-------------------------------------------------------------------------%
% Dynamic pole 3MUT
%-------------------------------------------------------------------------%
% col=hsv(length(t));
% for timeIdx = 3:1:length(t)
%     
%     % Define the time vector
%         tTemp = t(1:timeIdx);
% 
%         % General exponential elongation 
%         % We are looking to only one pole --> elongation / 2
%         l = @(T) ((exp(Kl*T)) + lc)/2;
% 
%         % General linear elongation 
%         %l = @(T) Kl*T + lc;
% 
%         % Septum Contribution to the Elongation (Feingold model - h(t))
%         lF = @(T) (dMax-dMax*( (tg-T)/(tg -tc)));
% 
%         % Septum constriction (Feingold model and 2D area model)
%         cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);
%         c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)));
% 
%         % Derivative(rate) elongation and constriction (Feingold model and 2D area model)
%         [dldt, d2ldt2] = dfdxc(tTemp, l(tTemp));
%         [dldtF, d2ldt2] = dfdxc(tTemp, lF(tTemp));
%         [dcdtF, d2cdt2] = dfdxc(tTemp, cF(tTemp));
%         [dcdt, d2cdt2] = dfdxc(tTemp, c(tTemp));
% 
%         % Choose an initial point
%         xd4 = 1;
%         yd4 = 2;
%         xdF4 = 1;
%         ydF4 = 2;
% 
%         % Build up the pole contour adding a costriction vector and an elongation
%         % vector depending on the PG insertion along the bacteria length
% 
%         % Initialize lTot
%         lTot = 0;
%         lTotF = 0;
%         for i = 2: length(tTemp) % cycle over the time steps
% 
%             % Cycle to grab the right elongation at each step 
%             % that has to be summed up to this new point
%             for dlIdx = i: length(tTemp) % from the elongation step i to the last
%                 % 2D model
%                 lTot = dldt(i-1)./i + lTot;
% 
%                 % Feingol model
%                 %lTotF = dldtF(i-1)./i + lTotF;
%                 %lTotF = (dldtF(i-1) + dldt(i-1))./i + lTotF;
%             end
% 
%             % Feingold model
%             xdF4(i) = xdF4(i-1) + dldtF(i-1)*(tTemp(i)-tTemp(i-1));
%             ydF4(i) = ydF4(i-1) + dcdtF(i-1);
%             % 2D model
%             xd4(i) = xd4(i-1) + lTot*(tTemp(i)-tTemp(i-1));
%             yd4(i) = yd4(i-1) + dcdt(i-1);
%         end % end cycle over the time step
% 
%         % contour top (left and right)
%         xd4TempTop = [xd4' - xd4(end) + x4(end); 2*xd4(end) - flipud(xd4') - xd4(end) + x4(end)];
%         yd4TempTop = [yd4'; flipud(yd4')];
%         % contour bottom (left and right)
%         xd4TempBot = [xd4' - xd4(end) + x4(end); 2*xd4(end) - flipud(xd4') - xd4(end) + x4(end)];
%         yd4TempBot = [-2*dMax - yd4'; flipud(-2*dMax - yd4')];
% 
%         % contour top (left and right)
%         xd4TempTopF = [xdF4' - xdF4(end) + xF4(end); 2*xdF4(end) - flipud(xdF4') + xF4(end)];
%         yd4TempTopF = [ydF4'; flipud(ydF4')];
%         % contour bottom (left and right)
%         xd4TempBotF = [xdF4' - xdF4(end) + xF4(end); 2*xdF4(end) - flipud(xdF4') + xF4(end)];
%         yd4TempBotF = [-2*dMax - ydF4'; flipud(-2*dMax - ydF4')];
% 
%         % contour with both right and left pole 
%         contourMUTAlld = [xd4TempTop, yd4TempTop, xd4TempBot, yd4TempBot];
%         contourMUTFAlld = [xd4TempTopF, yd4TempTopF, xd4TempBotF, yd4TempBotF];
% 
%         % contour of only the left pole
%         contourMUTHalfd = [xd4' - xd4(end) + x4(end), yd4', xd4' - xd4(end) + x4(end),- 2*dMax - yd4'];
%         contourMUTFHalfd = [xdF4'- xdF4(end) + xF4(end), ydF4', xdF4' - xdF4(end) + xF4(end),- 2*dMax - ydF4'];
% 
%         % DEBUG
%         subplot(1, 2, 2)
%         plot(contourMUTAlld(:, 1), contourMUTAlld(:, 2), '*', 'Color',col(i,:), 'markersize', 1),
%         hold on,
%         plot(contourMUTAlld(:, 3), contourMUTAlld(:, 4), '*', 'Color',col(i,:), 'markersize', 1),
%         axis equal
%         xlabel('length [nm]'),
%         ylabel('length [nm]'),
%         legend('3 MUT simulation'),
%         drawnow
% end % end dynamic pole
%-------------------------------------------------------------------------%

% figure, plot(x4, y4, '*');
% legend('bacteria contour')
% figure, plot(t, dcdt,'*')
% legend ('constrictio derivative')


lTot = 0;
lTotF = 0;
%-------------------------------------------------------------------------%
% MUT FtsW
%-------------------------------------------------------------------------%

%Kl = 0.014;
tc = 148;
tg = 360;
lc = 2180; % [nm]

t2 = (tc:1:tg);

l = @(T) ((exp(Kl*T)) + lc)/2;
lF = @(T)  (dMax-dMax*( (tg - T)/(tg - tc)));


c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)) );
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2);

% DEBUG
% figure, plot(t2, c(t2)),
% hold on
% plot(t2, l(t2))
% legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t2, l(t2));
[dldtF, d2ldt2] = dfdxc(t2, lF(t2));
[dcdt, d2cdt2] = dfdxc(t2, c(t2));
[dcdtF, d2cdt2] = dfdxc(t2, cF(t2));

x2=1;
y2=2;
xF2=1;
yF2=2;
for i = 2: length(t2)

    for dlIdx = i: length(t2) % from the elongation step i to the last
        % 2D model
        lTot = dldt(i-1)./i + lTot;
        % Feingol model
        %lTotF = dldtF(i-1)./i + lTotF;
    end
    x2(i) = x2(i-1) + lTot*(t2(i)-t2(i-1));
    y2(i) = y2(i-1) + dcdt(i-1);
    
    xF2(i) = xF2(i-1) + dldtF(i-1)*(t2(i)-t2(i-1));
    yF2(i) = yF2(i-1) + dcdtF(i-1);

end

% DEBUG
% figure, plot(t2, dcdt,'*')
% legend ('constriction derivative')
% 
% figure, plot(x1, y1, 'b*');
% hold on
% plot(x2, y2, 'r*');
% legend('bacteria contour tc=67 and tg=400', 'bacteria contour tc=93 and tg=290')


lTot = 0;
lTotF = 0;
%-------------------------------------------------------------------------%
% FtsI
%-------------------------------------------------------------------------%

% Kl = 0.014;
tc = 137;
tg = 337;
lc = 2343; % [nm]

t3 = (tc:1:tg);

l = @(T) ( exp(Kl*T) + lc )/2;
lF = @(T) (dMax-dMax*( (tg-T)/(tg -tc)));

c = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)) );
cF = @(T) dMax*wMax*sqrt( 1 - ((T-tc)/(tg-tc)).^2 );

% DEBUG
% figure, plot(t3, c(t3)),
% hold on
% plot(t3, l(t3))
% legend('constriction vs time', 'length vs time')


[dldt, d2ldt2] = dfdxc(t3, l(t3));
[dldtF, d2ldt2] = dfdxc(t3, lF(t3));
[dcdt, d2cdt2] = dfdxc(t3, c(t3));
[dcdtF, d2cdt2] = dfdxc(t3, cF(t3));

x3=1;
y3=2;
xF3=1;
yF3=2;

for i =2: length(t3)
    for dlIdx = i: length(t3) % from the elongation step i to the last
         % 2D model
        lTot = dldt(i-1)./i + lTot;
        % Feingol model
        %lTotF = dldtF(i-1)./i + lTotF;  
    end
    x3(i) = x3(i-1) + lTot*(t3(i)-t3(i-1));
    y3(i) = y3(i-1) + dcdt(i-1);
    xF3(i) = xF3(i-1) + dldtF(i-1)*(t3(i)-t3(i-1));
    yF3(i) = yF3(i-1) + dcdtF(i-1);
end

% DEBUG
% figure, plot(t3, dcdt,'*')
% legend ('constrictio derivative')

% 2D area - PG insertion model plot
figure, plot(x1, y1, 'b*'),
hold on
plot(x1, 2*y1(end)- y1, 'b*'),
hold all
plot(x2, y2 - (y1(end) - y2(end)), 'r*');
plot(x2, 2*y2(end)- y2  , 'r*'),
plot(x3, y3, 'g*');
plot(x3, 2*y3(end)- y3, 'g*'),
plot(x4, y4, 'k*');
plot(x4, 2*y4(end)- y4, 'k*'),
xlabel('[nm]')
ylabel('[nm]')
axis equal
title('2D area model - PG insertion')
legend('bacteria contour tc=67 and tg=400, WT', 'bacteria contour tc=67 and tg=400, WT',...
    'bacteria contour tc=148 and tg=360, FtsW(A246T)', 'bacteria contour tc=148 and tg=360, FtsW(A246T)',...
'bacteria contour tc=137 and tg=337 (FtsI I45V)', 'bacteria contour tc=137 and tg=337 (FtsI I45V)',...
'bacteria contour tc=118 and tg=269 3 MUT FtsW39', 'bacteria contour tc=95 and tg=253 FtsW39')

% 2D area - PG insertion model all septum plot
figure, plot(contourWTAll(:, 1), contourWTAll(:, 2), 'b*'),
hold on
plot(contourWTAll(:, 3), contourWTAll(:, 4), 'b*'),

plot(contourMUTWIAll(:, 1), contourMUTWIAll(:, 2), 'k*');
plot(contourMUTWIAll(:, 3), contourMUTWIAll(:, 4), 'k*'),

xlabel('[nm]')
ylabel('[nm]')
axis equal
title('2D area model - PG insertion')
legend('bacteria contour tc=67 and tg=400, WT', 'bacteria contour tc=67 and tg=400, WT',...
'bacteria contour tc=118 and tg=269 3 MUT FtsW39', 'bacteria contour tc=95 and tg=253 FtsW39')

% Faingold model plot
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
title('Feingold Model - PG Insertion')
legend('bacteria contour tc=67 and tg=400, WT', 'bacteria contour tc=67 and tg=400, WT',...
    'bacteria contour tc=148 and tg=360, FtsW(A246T)', 'bacteria contour tc=148 and tg=360, FtsW(A246T)',...
'bacteria contour tc=137 and tg=337 (FtsI I45V)', 'bacteria contour tc=137 and tg=337 (FtsI I45V)',...
'bacteria contour tc=118 and tg=269 FtsW39', 'bacteria contour tc=95 and tg=253 FtsW39')

% Save contours 2D model in 2 coloumns
wtContour = [[x1 x1]; [y1 2*y1(end)- y1]]';
mutFtsWContour = [[x2 x2]; [y2 + (y1(end) - y2(end)) 2*y2(end)- y2 + (y1(end) - y2(end))]]';
mutFtsIContour = [[x3 x3]; [y3 + (y1(end) - y3(end)) 2*y3(end)- y3 + (y1(end) - y3(end))]]';
mutFtsWIContour = [[x4 x4]; [y4 + (y1(end) - y4(end)) 2*y4(end)- y4 + (y1(end) - y4(end))]]';

% Save contours Feingold model in 2 coloumns
wtContourF = [[xF1 xF1]; [yF1 2*yF1(end)- yF1]]';
mutFtsWContourF = [[xF2 xF2]; [yF2 + (yF1(end) - yF2(end)) 2*yF2(end)- yF2 + (yF1(end) - yF2(end))]]';
mutFtsIContourF = [[xF3 xF3]; [yF3 + (yF1(end) - yF3(end)) 2*yF3(end)- yF3 + (yF1(end) - yF3(end))]]';
mutFtsWIContourF = [[xF4 xF4]; [yF4 + (yF1(end) - yF4(end)) 2*yF4(end)- yF4 + (yF1(end) - yF4(end))]]';

save('G:\Anna\SimuletorResults\MatlabVariables\2DAreaModContour.mat', 'wtContour', 'mutFtsWContour', 'mutFtsIContour', 'mutFtsWIContour');
save('G:\Anna\SimuletorResults\MatlabVariables\FeingModContour.mat', 'wtContourF', 'mutFtsWContourF', 'mutFtsIContourF', 'mutFtsWIContourF');

% Save contours 2D model in 4 coloumns
save('G:\Anna\SimuletorResults\MatlabVariables\2DAreaModContour4Column.mat', 'contourWTHalf', 'contourMUTWIHalf');
save('G:\Anna\SimuletorResults\MatlabVariables\FeingModContour4Column.mat', 'contourWTFHalf', 'contourMUTWIFHalf');

% Debug
% figure, plot(wtContour(:,1), wtContour(:,2),'*'), hold on,
% plot(mutFtsWContour(:,1), mutFtsWContour(:,2),'*'),
% plot(mutFtsIContour(:,1), mutFtsIContour(:,2),'*'),
% plot(mutFtsWIContour(:,1), mutFtsWIContour(:,2),'*'),

% Debug
% figure, plot(wtContourF(:,1), wtContourF(:,2),'*'), hold on,
% plot(mutFtsWContourF(:,1), mutFtsWContourF(:,2),'*'),
% plot(mutFtsIContourF(:,1), mutFtsIContourF(:,2),'*'),
% plot(mutFtsWIContourF(:,1), mutFtsWIContourF(:,2),'*'),

