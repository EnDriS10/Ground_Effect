% --- Gráfica de Datos de Sustentación vs Altura (Con Tendencia) ---
clear; clc; close all;

% 1. Datos originales de la tabla
h_mm = [1.74, 11.74, 21.70, 94.6, 171.74];
L_N = [16.5106, 12.4192, 10.7362, 5.6709, 5.4512];

% 2. Crear un vector denso para la curva de tendencia suave
h_tendencia = linspace(min(h_mm), max(h_mm), 200);

% Calculamos la línea de tendencia suave (pchip evita oscilaciones irreales)
L_tendencia = pchip(h_mm, L_N, h_tendencia);

% 3. Crear la gráfica 2D Estilo Paper
figure('Color', 'w', 'Position', [200, 200, 700, 500]);
hold on; % Permite dibujar varias cosas en la misma gráfica

% Dibujamos la línea de tendencia suave en azul puro ('b'), igual que la otra gráfica
plot(h_tendencia, L_tendencia, 'b-', 'LineWidth', 2.5);

% Dibujamos los puntos de datos sueltos encima (círculos negros para que destaquen)
plot(h_mm, L_N, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'LineWidth', 1);

% 4. Etiquetas y título
title('Efecto Suelo: Sustentación vs. Altura', 'FontSize', 14,'Color','k');
xlabel('Altura sobre el suelo, h (mm)', 'FontSize', 12, 'Color', 'k');
ylabel('Sustentación, L (N)', 'FontSize', 12, 'Color', 'k');

% 5. Configuración del fondo y los ejes
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'FontSize', 11, 'GridAlpha', 0.15, 'LineWidth', 1);
grid on;
box on;

% 6. Leyenda para que quede súper profesional en el paper
lgd = legend('Tendencia Suave', 'Datos Medidos', 'Location', 'northeast', 'FontSize', 11);

% Personalización de la leyenda (Fondo blanco, texto y borde negros)
lgd.Color = 'w';       % 'w' es white (blanco) para el fondo
lgd.TextColor = 'k';   % 'k' es black (negro) para las letras
lgd.EdgeColor = 'k';   % Opcional: Borde de la cajita en negro para que resalte
% Ajustar los límites de los ejes
xlim([0, 180]);
ylim([4, 18]);
hold off;

% exportgraphics(gcf, 'sustentacion_tendencia.pdf', 'ContentType', 'vector');