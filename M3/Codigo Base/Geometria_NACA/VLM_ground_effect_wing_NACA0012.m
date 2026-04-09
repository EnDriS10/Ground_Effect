clear; clc; close all;

%% Parámetros
b  = 10.0;
c  = 1.0;
N_x = 10;
N_y = 20;   % menos para visualizar mejor
L_wake = 5*c;

%% -------------------------------
% MALLA CON DISTRIBUCIÓN COSENO
% -------------------------------

% Cuerda
beta_c   = linspace(0, pi, N_x + 1);
x_edges  = c/2 * (1 - cos(beta_c));
dx_panel = diff(x_edges);

x_vort = x_edges(1:N_x) + dx_panel/4;
x_ctrl = x_edges(1:N_x) + 3*dx_panel/4;

% Envergadura
beta_s  = linspace(0, pi, N_y + 1);
y_edges = -b/2 * cos(beta_s);
y_ctrl  = 0.5 * (y_edges(1:N_y) + y_edges(2:end));

%% -------------------------------
% GEOMETRÍA NACA 0012
% -------------------------------
% Función que devuelve z = ±t(x) sobre la cuerda normalizada x/c ∈ [0,1]
% Ecuación NACA 4 dígitos: t(x/c) = 5*0.12*c * (a0*sqrt(xi)-a1*xi-a2*xi²+a3*xi³-a4*xi⁴)
% con xi = x/c y coeficientes estándar NACA.

naca0012 = @(x_abs) naca_thickness(x_abs, c, 0.12);

%% -------------------------------
% FIGURA
% -------------------------------
figure('Color','w'); hold on; grid on;
view(3)

%% 1) Dibujar superficie del ala — NACA 0012
for j = 1:N_y
    for i = 1:N_x
        x1 = x_edges(i);
        x2 = x_edges(i+1);
        y1 = y_edges(j);
        y2 = y_edges(j+1);

        % Ordenadas z del extradós en los 4 vértices del panel
        z11 = naca0012(x1);   % (x1, y1) extradós
        z21 = naca0012(x2);   % (x2, y1)
        z22 = naca0012(x2);   % (x2, y2) — perfil constante en envergadura
        z12 = naca0012(x1);   % (x1, y2)

        % Panel del extradós (cara superior)
        fill3([x1 x2 x2 x1], [y1 y1 y2 y2], [z11 z21 z22 z12], ...
              [0.75 0.85 0.95], 'EdgeColor','k', 'FaceAlpha', 0.9);

        % Panel del intradós (cara inferior, z negativa)
        fill3([x1 x2 x2 x1], [y1 y1 y2 y2], [-z11 -z21 -z22 -z12], ...
              [0.65 0.75 0.85], 'EdgeColor','k', 'FaceAlpha', 0.7);
    end
end

%% 2) Dibujar vórtices ligados (líneas azules) — en el plano medio z=0
for j = 1:N_y
    for i = 1:N_x
        xv = x_vort(i);
        y1 = y_edges(j);
        y2 = y_edges(j+1);
        zv = naca0012(xv);   % z del extradós en xv

        % Vórtice ligado dibujado en la línea media (z ≈ 0 en VLM)
        plot3([xv xv], [y1 y2], [0 0], 'b', 'LineWidth', 2);
    end
end

%% 3) Dibujar puntos de control (rojo) — en el plano z=0 (camber line plana)
for j = 1:N_y
    for i = 1:N_x
        xc = x_ctrl(i);
        yc = y_ctrl(j);

        plot3(xc, yc, 0, 'ro', 'MarkerSize', 4, 'MarkerFaceColor','r');
    end
end

%% 4) Dibujar estela (wake) — sale del TE (x = c) a z = 0
for j = 1:N_y
    for i = 1:N_x
        xv = x_vort(i);
        y1 = y_edges(j);
        y2 = y_edges(j+1);

        % líneas hacia atrás desde cada vórtice
        plot3([xv xv+L_wake], [y1 y1], [0 0], 'k--');
        plot3([xv xv+L_wake], [y2 y2], [0 0], 'k--');
    end
end

%% 5) Dibujar perfil NACA 0012 en la sección central (y=0) y en los tips
xi_plot  = linspace(0, c, 200);
z_extrad = arrayfun(naca0012, xi_plot);

% Sección central
fill3(xi_plot, zeros(size(xi_plot)), z_extrad, 'b', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.1 0.3 0.7], 'LineWidth', 1.8);
fill3(xi_plot, zeros(size(xi_plot)), -z_extrad, 'b', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.1 0.3 0.7], 'LineWidth', 1.8);

% Tip izquierdo
fill3(xi_plot, -b/2*ones(size(xi_plot)), z_extrad, 'r', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.7 0.1 0.1], 'LineWidth', 1.8);
fill3(xi_plot, -b/2*ones(size(xi_plot)), -z_extrad, 'r', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.7 0.1 0.1], 'LineWidth', 1.8);

% Tip derecho
fill3(xi_plot, b/2*ones(size(xi_plot)), z_extrad, 'r', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.7 0.1 0.1], 'LineWidth', 1.8);
fill3(xi_plot, b/2*ones(size(xi_plot)), -z_extrad, 'r', 'FaceAlpha', 0.0, ...
      'EdgeColor',[0.7 0.1 0.1], 'LineWidth', 1.8);

%% 6) Estética
xlabel('x (cuerda)');
ylabel('y (envergadura)');
zlabel('z (espesor)');
title('Geometría VLM — Ala NACA 0012: paneles, vórtices y wake');

axis equal
xlim([0 c + L_wake])
ylim([-b/2 b/2])
zlim([-0.15*c  0.15*c])

legend({'Extradós','Intradós','Vórtices ligados','Puntos de control','Wake'}, ...
        'Location','best')

% =========================================================================
%  FUNCIÓN LOCAL: espesor NACA 4 dígitos simétrico
% =========================================================================
function z = naca_thickness(x_abs, c, t_ratio)
% Devuelve la ordenada z del extradós (z > 0) del perfil NACA 4 dígitos
% simétrico para una posición x_abs (en metros, 0 ≤ x_abs ≤ c).
%
%   t_ratio : fracción de espesor máximo (ej. 0.12 para NACA 0012)
%   Coeficientes NACA estándar (cierre exacto en TE con a4 modificado):

    xi   = x_abs / c;                  % coordenada normalizada [0, 1]
    xi   = max(0, min(1, xi));          % clamp por seguridad

    a0 =  0.2969;
    a1 = -0.1260;
    a2 = -0.3516;
    a3 =  0.2843;
    a4 = -0.1015;   % −0.1036 para TE abierto; −0.1015 para TE cerrado

    z = 5 * t_ratio * c * ...
        (a0*sqrt(xi) + a1*xi + a2*xi^2 + a3*xi^3 + a4*xi^4);
end
