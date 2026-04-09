close all; clear; clc;

b = 10;           % Envergadura (m)
Gamma = 25;       % Circulación (m^2/s)
V_viento = 15;       % Velocidad del aire (m/s)
L_estela = 80;    % Largo de la estela (m)
rho = 1.225;      % Densidad del aire (kg/m^3)
h = linspace(1,20,100);
L_numerico = zeros(size(h));
for i = 1:length(h)
    L_numerico(i)=calcular_sustentacion(h(i),b,Gamma,V_viento,rho,L_estela);
end
% 3. Gráfica 2D Estilo Paper
figure('Color', 'w', 'Position', [200, 200, 700, 500]);
plot(h, L_numerico, 'b-', 'LineWidth', 2.5); % Línea azul, continua y gruesa

% Etiquetas y título
title('Efecto Suelo: Sustentación vs. Altura', 'FontSize', 14,'Color','k');
xlabel('Altura sobre el suelo, h (m)', 'FontSize', 12, 'Color', 'k');
ylabel('Sustentación Numérica, L (N)', 'FontSize', 12, 'Color', 'k');

% Configuración del fondo y los ejes
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'FontSize', 11, 'GridAlpha', 0.15);
grid on;
box on; % Pone un recuadro cerrado alrededor de la gráfica


function [L_numerico] = calcular_sustentacion(h, b, Gamma, V_viento, rho, L_estela)
    % calcular_sustentacion: Calcula la fuerza de sustentación (Lift) de un ala 
    % en efecto suelo, devolviendo tanto el modelo teórico como el numérico.
    
    % --- 1. CÁLCULO ANALÍTICO (Teórico de Prandtl) ---
    w_inf = Gamma / (-pi * b); 
    L_inf = rho * V_viento * Gamma * b;
    
    % Factor empírico del efecto suelo
    factor_h = (16 * (h/b)^2) / (1 + 16 * (h/b)^2);
    w_total_teorico = w_inf * factor_h;
    
    % Sustentación Teórica Final
    L_teorico = L_inf * (w_inf / w_total_teorico);

    % --- 2. CÁLCULO NUMÉRICO (Biot-Savart + Trefftz) ---
    % Definir filamentos (Real y Espejo)
    real_fil = {[L_estela, -b/2, h; 0, -b/2, h], ...
                [0, -b/2, h; 0, b/2, h], ...
                [0, b/2, h; L_estela, b/2, h]};
    espj_fil = {[L_estela, -b/2, -h; 0, -b/2, -h], ...
                [0, -b/2, -h; 0, b/2, -h], ...
                [0, b/2, -h; L_estela, b/2, -h]};
                
    filamentos = [real_fil, espj_fil];
    gammas = [Gamma, Gamma, Gamma, -Gamma, -Gamma, -Gamma];
    
    % Punto de evaluación en la estela lejana
    P = [40, 0, h];
    W_punto = 0; 
    
    % Evaluación analítica de Biot-Savart en el punto P
    for f = 1:6
        seg = filamentos{f}; 
        G = gammas(f);
        A = seg(1,:); 
        B = seg(2,:);
        
        r1 = P - A; 
        r2 = P - B; 
        r0 = B - A;
        
        c = cross(r1, r2); 
        mag_c2 = sum(c.^2);
        
        if mag_c2 > 1e-6 
            V_ind = (G / (4 * pi)) * (c / mag_c2) * (dot(r1, r0)/norm(r1) - dot(r2, r0)/norm(r2));
            W_punto = W_punto + V_ind(3);
        end
    end
    
    % Extracción de datos numéricos (Teorema de Trefftz)
    w_num = W_punto / 2;
    L_numerico = L_inf * (w_inf / w_num);
end