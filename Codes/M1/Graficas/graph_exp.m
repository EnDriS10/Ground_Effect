% --- Gráfica de Datos Experimentales (Sin el primer dato) ---
clear; clc; close all;

% 1. Datos extraídos (Dato de h=25 eliminado)
h_mm  = [45, 55, 70, 80, 115];
err_h = [1, 1, 1, 1, 1]; % Incertidumbre h (mm)

L_mN  = [520, 410, 390, 330, 300];
err_L = [40, 40, 40, 30, 40]; % Incertidumbre L (mN)

% 2. CÁLCULO DE LA LÍNEA DE TENDENCIA SUAVE
% Usamos una regresión polinómica de 2º orden para 5 puntos
p = polyfit(h_mm, L_mN, 2); 

% Crear un vector denso para dibujar la curva suave
h_tendencia = linspace(min(h_mm), max(h_mm), 200);
L_tendencia = polyval(p, h_tendencia); 

% 3. Crear la gráfica 2D Estilo Paper
figure('Color', 'w', 'Position', [200, 200, 700, 500]);
hold on;

% Dibujamos la línea de tendencia suave en azul puro ('b')
plot(h_tendencia, L_tendencia, 'b-', 'LineWidth', 2.5);

% Dibujamos los puntos experimentales CON barras de error ENCIMA
errorbar(h_mm, L_mN, err_L, err_L, err_h, err_h, 'ko', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'CapSize', 4);

% 4. Etiquetas y título
title('Efecto Suelo Experimental: Sustentación vs. Altura', 'FontSize', 14, 'Color','k');
xlabel('Altura sobre el suelo, h (mm)', 'FontSize', 12, 'Color', 'k');
ylabel('Fuerza de sustentación, L (mN)', 'FontSize', 12, 'Color', 'k');

% 5. Configuración del fondo y los ejes
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'FontSize', 11, 'GridAlpha', 0.15, 'LineWidth', 1);
grid on;
box on;

% 6. Leyenda personalizada (Fondo blanco, texto