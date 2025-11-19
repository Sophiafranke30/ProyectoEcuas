clear; clc; close all;

%% CONFIGURACIÓN GENERAL
h = 0.1; % Tamaño de paso (se puede variar para el estudio de convergencia |valores recomendados -> h = 0.05 y h = 0.01)
tspan = [0, 10]; % Intervalo de tiempo de 0 a 10

%% CASO 1: ED DE PRIMER ORDEN
% Ecuación: dy/dx = 3x - y, y(0) = 2
% Solución Analítica (resolviendo factor integrante): y = 3x - 3 + 5*exp(-x)

disp('--- Resolviendo ED Primer Orden ---');
f1 = @(x, y) 3*x - y;
y0_1 = 2;

% Solución Numérica
[t1_eu, y1_eu] = metodo_euler(f1, tspan, y0_1, h);
[t1_rk, y1_rk] = metodo_rk4(f1, tspan, y0_1, h);

% Solución Analítica
y_exact1 = 3*t1_rk - 3 + 5*exp(-t1_rk);

% Gráfica
figure(1);
plot(t1_eu, y1_eu, 'r--', 'LineWidth', 1.5); hold on;
plot(t1_rk, y1_rk, 'b-.', 'LineWidth', 1.5);
plot(t1_rk, y_exact1, 'k', 'LineWidth', 1);
legend('Euler', 'RK4', 'Analítica');
title('ED Primer Orden: dy/dx = 3x - y');
xlabel('x'); ylabel('y'); grid on;

%% CASO 2: ED DE SEGUNDO ORDEN
% Ecuación: y'' + 4y' + 3y = 0, y(0)=0, y'(0)=1
% Conversión a sistema: u1 = y, u2 = y'
% u1' = u2
% u2' = -3*u1 - 4*u2
% Solución Analítica: y = 0.5*exp(-x) - 0.5*exp(-3*x)

disp('--- Resolviendo ED Segundo Orden ---');
f2 = @(t, u) [u(2); -3*u(1) - 4*u(2)];
y0_2 = [0; 1]; % [y(0); y'(0)]

[t2_eu, y2_eu] = metodo_euler(f2, tspan, y0_2, h);
[t2_rk, y2_rk] = metodo_rk4(f2, tspan, y0_2, h);

% Solución Analítica para y (que es u1)
y_exact2 = 0.5*exp(-t2_rk) - 0.5*exp(-3*t2_rk);

% Gráfica
figure(2);
plot(t2_eu, y2_eu(:,1), 'r--', 'LineWidth', 1.5); hold on; % Solo graficamos y (columna 1)
plot(t2_rk, y2_rk(:,1), 'b-.', 'LineWidth', 1.5);
plot(t2_rk, y_exact2, 'k', 'LineWidth', 1);
legend('Euler (y)', 'RK4 (y)', 'Analítica');
title('ED Segundo Orden: y'''' + 4y'' + 3y = 0');
xlabel('t'); ylabel('y(t)'); grid on;

%% CASO 3: SISTEMA LINEAL 2x2
% dx/dt = 2x + y
% dy/dt = x + y
% x(0)=1, y(0)=0

disp('--- Resolviendo Sistema Lineal 2x2 ---');
f3 = @(t, v) [2*v(1) + v(2); v(1) + v(2)]; % v(1) es x, v(2) es y
y0_3 = [1; 0];

[t3_eu, y3_eu] = metodo_euler(f3, tspan, y0_3, h);
[t3_rk, y3_rk] = metodo_rk4(f3, tspan, y0_3, h);

% Gráfica (x vs t y y vs t)
figure(3);
subplot(2,1,1);
plot(t3_eu, y3_eu(:,1), 'r--', t3_rk, y3_rk(:,1), 'b-.');
legend('Euler x(t)', 'RK4 x(t)'); title('Variable x(t)'); grid on;
subplot(2,1,2);
plot(t3_eu, y3_eu(:,2), 'r--', t3_rk, y3_rk(:,2), 'b-.');
legend('Euler y(t)', 'RK4 y(t)'); title('Variable y(t)'); grid on;

%% CASO 4: SISTEMA NO LINEAL (Van der Pol)
% x' = y
% y' = mu*(1-x^2)*y - x
% x(0)=2, y(0)=0, mu > 0 (Usaremos mu=1)

disp('--- Resolviendo Van der Pol (No Lineal) ---');
mu = 1; 
f4 = @(t, z) [z(2); mu*(1 - z(1)^2)*z(2) - z(1)];
y0_4 = [2; 0];
tspan_vdp = [0, 20]; % Necesitamos más tiempo para ver la oscilación

[t4_eu, y4_eu] = metodo_euler(f4, tspan_vdp, y0_4, h);
[t4_rk, y4_rk] = metodo_rk4(f4, tspan_vdp, y0_4, h);

% Gráfica de Plano de Fase (y vs x)
figure(4);
plot(y4_eu(:,1), y4_eu(:,2), 'r--'); hold on;
plot(y4_rk(:,1), y4_rk(:,2), 'b-');
legend('Euler (Plano Fase)', 'RK4 (Plano Fase)');
title(['Oscilador de Van der Pol (\mu = ' num2str(mu) ')']);
xlabel('x'); ylabel('y (o dx/dt)'); grid on;
