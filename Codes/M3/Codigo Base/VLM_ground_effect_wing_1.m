% =========================================================================
%  VLM - Visualización con Ángulo de Ataque (alpha) y Efecto Suelo
% =========================================================================
clear; clc; close all;

%% 1 — PARÁMETROS
b        = 10.0;        
c        = 1.0;        
alpha    = deg2rad(5); % Ángulo de ataque del código original
N_x      = 6;          
N_y      = 14;         
L_wake   = 3 * c;      
h_rel    = 0.4;        
h_pivot  = h_rel * c;  % Altura en el borde de ataque (pivot)

%% 2 — GENERACIÓN DE MALLA
[Xv, Xc, Yc, YL, YR, Dx, x_edges, y_edges, x_ctrl, y_ctrl, dx_panel] = Malla_VLM(b, c, N_x, N_y);

%% 3 — CONFIGURACIÓN DE FIGURA
figure('Color', 'w', 'Name', 'VLM Geometry with Alpha');
hold on; grid on; box on;
view(35, 25); axis equal;
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

% Colores
col_real = [0.2, 0.6, 0.9]; 
col_img  = [0.8, 0.2, 0.2]; 
col_vort = [0.0, 0.0, 0.8];

%% 4 — TRANSFORMACIÓN Y DIBUJO
% Función anónima para calcular Z con rotación: z = h_pivot - x * sin(alpha)
calc_z = @(x, h_val, a_val) h_val - x * sin(a_val);

% --- Dibujo de Paneles (Ala Real e Imagen) ---
for j = 1:N_y
    for i = 1:N_x
        px = [x_edges(i), x_edges(i+1), x_edges(i+1), x_edges(i)];
        py = [y_edges(j), y_edges(j), y_edges(j+1), y_edges(j+1)];
        
        % Coordenadas Z rotadas
        pz_real = calc_z(px, h_pivot, alpha);
        pz_img  = -pz_real; % Espejo exacto respecto a z=0
        
        % Patch Real
        patch(px, py, pz_real, col_real, 'FaceAlpha', 0.3, 'EdgeColor', [0.3 0.3 0.3]);
        % Patch Imagen
        patch(px, py, pz_img, col_img, 'FaceAlpha', 0.1, 'EdgeColor', [0.7 0.7 0.7], 'LineStyle', ':');
    end
end

% --- Vórtices, Puntos de Control y Estela ---
for k = 1:length(Xv)
    % Z para el vórtice y punto de control
    zv = calc_z(Xv(k), h_pivot, alpha);
    zc = calc_z(Xc(k), h_pivot, alpha);
    
    % 1. Vórtice Ligado (Bound Vortex)
    plot3([Xv(k) Xv(k)], [YL(k) YR(k)], [zv zv], 'Color', col_vort, 'LineWidth', 2);
    
    % 2. Punto de Control (Punto Rojo)
    plot3(Xc(k), Yc(k), zc, 'r.', 'MarkerSize', 10);
    
    % 3. Estela (Wake) - Sale desde los extremos del vórtice ligado
    % Nota: La estela se alinea con la corriente libre (eje X)
    plot3([Xv(k) Xv(k)+L_wake], [YL(k) YL(k)], [zv zv], 'Color', [0.4 0.4 0.4 0.5]);
    plot3([Xv(k) Xv(k)+L_wake], [YR(k) YR(k)], [zv zv], 'Color', [0.4 0.4 0.4 0.5]);
end

%% 5 — ELEMENTOS DE REFERENCIA
% Plano del Suelo (z=0)
fill3([-0.5, c+L_wake, c+L_wake, -0.5], [-b, -b, b, b], [0 0 0 0], ...
      [0.95 0.95 0.95], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

% Etiquetas y Títulos
title(sprintf('Geometría VLM : \\alpha = %.1f°, h/c = %.2f', rad2deg(alpha), h_rel));
xlabel('x (Cuerda)'); ylabel('y (Envergadura)'); zlabel('z (Altura)');

% Leyenda
dummy_real = plot(nan, nan, 'Color', col_real, 'LineWidth', 5);
dummy_vort = plot(nan, nan, 'Color', col_vort, 'LineWidth', 2);
legend([dummy_real, dummy_vort], {'Ala (Incidencia \alpha)', 'Vórtices Ligados'}, 'Location', 'northeast');

% Ajustar límites
zlim([0, 2*h_rel]);
%%
%  =========================================================================
%  FUNCIONES LOCALES
% =========================================================================
function [Xv, Xc, Yc, YL, YR, Dx, x_edges, y_edges, x_ctrl, y_ctrl, dx_panel] = Malla_VLM(b, c, N_x, N_y)
    N_tot = N_x * N_y;
    beta_c = linspace(0, pi, N_x + 1);
    x_edges = c/2 * (1 - cos(beta_c));
    dx_panel = diff(x_edges);
    x_vort = x_edges(1:N_x) + dx_panel / 4;
    x_ctrl = x_edges(1:N_x) + 3*dx_panel / 4;
    beta_s = linspace(0, pi, N_y + 1);
    y_edges = -b/2 * cos(beta_s);
    y_ctrl = 0.5 * (y_edges(1:N_y) + y_edges(2:N_y+1));
    
    Xv = zeros(N_tot, 1); Xc = zeros(N_tot, 1); Yc = zeros(N_tot, 1);
    YL = zeros(N_tot, 1); YR = zeros(N_tot, 1); Dx = zeros(N_tot, 1);
    for j = 1:N_y
        for i = 1:N_x
            k = (j-1)*N_x + i;
            Xv(k) = x_vort(i); Xc(k) = x_ctrl(i); Yc(k) = y_ctrl(j);
            YL(k) = y_edges(j); YR(k) = y_edges(j+1); Dx(k) = dx_panel(i);
        end
    end
end