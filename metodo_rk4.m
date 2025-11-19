function [t, y] = metodo_rk4(f, tspan, y0, h)
    % Implementación de Runge-Kutta de 4to Orden
    
    t = tspan(1):h:tspan(2);
    N = length(t);
    m = length(y0);
    y = zeros(m, N);
    y(:, 1) = y0;
    
    for i = 1:N-1
        ti = t(i);
        yi = y(:, i);
        
        k1 = f(ti, yi);
        k2 = f(ti + 0.5*h, yi + 0.5*h*k1);
        k3 = f(ti + 0.5*h, yi + 0.5*h*k2);
        k4 = f(ti + h, yi + h*k3);
        
        y(:, i+1) = yi + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    y = y';
end
