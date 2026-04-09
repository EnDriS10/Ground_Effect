% --- Script: Graficar Error 3D f(h, b) ---
clear; clc; close all;

% 1. Definir los vectores
b_vec = linspace(1, 10, 100);
h_vec = linspace(0.1, 20, 100);
V_inf = 15; % Mantenemos V_inf constante

% 2. Crear las mallas 2D
[B, H] = meshgrid(b_vec, h_vec);

% Matriz vacía para guardar los resultados del error
ERROR_Z = zeros(size(B));

disp('Calculando 10,000 puntos de la superficie de error... espere un momento.');

% 3. Bucle para evaluar la función en cada punto (h, b)
for i = 1:size(B, 1)
    for j = 1:size(B, 2)
        b_actual = B(i,j);
        h_actual = H(i,j);
        
        % Llamamos a la función local M1
        ERROR_Z(i,j) = M1(h_actual, b_actual);
    end
end

disp('Cálculo terminado. Generando gráfico...');

% =========================================================================
% 4. Visualización 3D de la superficie (ESTILO PAPER)
% =========================================================================
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
surf(B, H, ERROR_Z, 'EdgeColor', 'none'); 

% Usamos parula, que es el estándar científico de MATLAB
colormap(parula); 
cb = colorbar;
cb.Color = 'k';
ylabel(cb, 'Error (%)', 'FontSize', 12, 'Color', 'k');

% Etiquetas (Eliminamos el título, eso va en LaTeX)
xlabel('Envergadura b (m)', 'FontSize', 12, 'Color', 'k');
ylabel('Altura h (m)', 'FontSize', 12, 'Color', 'k');
zlabel('Diferencia de Sust. (%)', 'FontSize', 12, 'Color', 'k');

% Forzamos fondo blanco, texto negro y malla visible pero sutil
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
         'GridColor', 'k', 'GridAlpha', 0.15, 'FontSize', 11);
         
view([-40, 35]); 
grid on;

% ¡NUEVO! Exportar automáticamente como PDF vectorial
exportgraphics(gcf, 'grafica_error.pdf', 'ContentType', 'vector');


% =========================================================================
% FUNCIONES LOCALES (Deben ir siempre al final del archivo)
% =========================================================================

function error_numerico = M1(h, b)
    % M1: Simula el efecto suelo aerodinámico.
    
    % --- ¡CORRECCIÓN CRÍTICA! ---
    % Control de gráficos a 'false' para que no abra 10,000 ventanas
    mostrar_graficos = false; 
    
    if mostrar_graficos
        figure('Color', 'k', 'Position', [100 100 1000 700]);
    end
    
    % --- 1. PARÁMETROS FÍSICOS ---
    V_inf = 15;
    Gamma = 25;       % Circulación (m^2/s)
    L_estela = 80;    % Largo de la estela (m)
    rho = 1.225;      % Densidad del aire (kg/m^3)
    
    % --- 2. CÁLCULO ANALÍTICO (Fórmulas Teóricas) ---
    w_inf = Gamma / (-pi * b); 
    
    % Factor empírico del efecto suelo
    factor_h = (16 * (h/b)^2) / (1 + 16 * (h/b)^2);
    w_total_teorico = w_inf * factor_h;
    
    L_inf = rho * V_inf * Gamma * b;
    L_total_teorico = L_inf * (w_inf / w_total_teorico);
    
    % --- 3. CÁLCULO NUMÉRICO Y ERROR ---
    L_numerico = calcular_L_numerico(h, b, Gamma, V_inf, rho, L_estela);
    error_numerico = 100 * abs((L_total_teorico - L_numerico) / L_total_teorico);
end

function L_num = calcular_L_numerico(h_actual, b, Gamma, V_inf, rho, L_estela)
    % 1. Parámetros base
    w_inf = Gamma / (-pi * b); 
    L_inf = rho * V_inf * Gamma * b;
    
    % 2. Definir filamentos (Real y Espejo)
    real_fil = {[L_estela, -b/2, h_actual; 0, -b/2, h_actual], ...
                [0, -b/2, h_actual; 0, b/2, h_actual], ...
                [0, b/2, h_actual; L_estela, b/2, h_actual]};
    espj_fil = {[L_estela, -b/2, -h_actual; 0, -b/2, -h_actual], ...
                [0, -b/2, -h_actual; 0, b/2, -h_actual], ...
                [0, b/2, -h_actual; L_estela, b/2, -h_actual]};
                
    filamentos = [real_fil, espj_fil];
    gammas = [Gamma, Gamma, Gamma, -Gamma, -Gamma, -Gamma];
    
    % 3. Punto de evaluación en la estela lejana
    P = [40, 0, h_actual];
    W_punto = 0; 
    
    % 4. Biot-Savart evaluado solo en el punto P (¡Esto lo habías borrado!)
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
    
    % 5. Extraer dato y aplicar Teorema de Trefftz
    w_num = W_punto / 2;
    
    % 6. Calcular sustentación final basada en el downwash numérico
    L_num = L_inf * (w_inf / w_num);
end