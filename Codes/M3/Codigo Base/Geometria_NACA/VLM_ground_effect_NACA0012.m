% =========================================================================
%  Vortex Lattice Method (VLM) con Efecto Suelo — Perfil NACA 0012
% =========================================================================

clear; clc; close all;

%%
% =========================================================================
%  BLOQUE 1 — PARÁMETROS
% =========================================================================

% Geometría
b        = 10.0;        % Envergadura [m]
c        = 1.0;        % Cuerda [m]

% Vuelo
V_inf    = 15;         % Velocidad de la corriente libre [m/s]
alpha    = deg2rad(5); % Ángulo de ataque [rad]
rho      = 1.225;      % Densidad del aire [kg/m³]
S_ref    = b * c;      % Área de referencia [m²]

[AR , e_oswald] = Oswald(b,c);
fprintf('Factor de Oswald calculado: e = %.3f (AR=%.1f)\n', e_oswald, AR);

% Discretización
N_x      = 10;
N_y      = 48;
N_tot    = N_x * N_y;

% Estela
L_wake   = 150 * c;

% Alturas sobre el suelo h/c
h_ratio  = [0.10, 0.15, 0.20, 0.25, 0.30,0.35,0.40,0.45, 0.50,0.55, ... 
    0.60,0.65,0.70,0.75,0.80,0.85,0.90, 1.00,1.25, 1.50,1.75, 2.00,2.50,3.00];
N_h      = length(h_ratio);

fprintf('VLM Efecto Suelo NACA 0012 | AR=%.1f | Nx=%d Ny=%d | Ntot=%d\n', ...
        AR, N_x, N_y, N_tot);
fprintf('Estela: %.0f×c | Oswald: e=%.3f\n\n', L_wake/c, e_oswald);

%%
% =========================================================================
%  BLOQUE 2 — MALLA VLM CON DISTRIBUCIÓN COSENO + GEOMETRÍA NACA 0012
% =========================================================================

% -- Cuerda: distribución coseno
beta_c   = linspace(0, pi, N_x + 1);
x_edges  = c/2 * (1 - cos(beta_c));          % [0 → c]
dx_panel = diff(x_edges);
x_vort   = x_edges(1:N_x) + dx_panel / 4;    % 1/4 cuerda local → vórtice
x_ctrl   = x_edges(1:N_x) + 3*dx_panel / 4;  % 3/4 cuerda local → control

% -- Envergadura: distribución coseno
beta_s   = linspace(0, pi, N_y + 1);
y_edges  = -b/2 * cos(beta_s);
y_ctrl   = 0.5 * (y_edges(1:N_y) + y_edges(2:N_y+1));

% -- Ordenadas z y pendiente dz/dx del perfil NACA 0012 en cada posición x
%    El VLM trata el ala como "camber surface"; para NACA 0012 (simétrico)
%    la camber line es z=0, pero el ángulo efectivo de la normal se obtiene
%    de la pendiente de la línea de espesor medio → dz/dx = 0 para perfil
%    simétrico. Sin embargo, desplazamos los puntos de control a la
%    superficie del extradós (z_ctrl > 0) para reflejar la geometría real.

z_vort = arrayfun(@(x) naca_z(x, c, 0.12), x_vort);   % z en línea 1/4 (extradós)
z_ctrl = arrayfun(@(x) naca_z(x, c, 0.12), x_ctrl);   % z en línea 3/4 (extradós)
dzdx   = arrayfun(@(x) naca_dzdx(x, c, 0.12), x_ctrl); % pendiente en ctrl

% Normal al extradós: n = (-dz/dx, 0, 1) normalizada en cada panel
% La condición de impermeabilidad es: V·n = 0
% RHS_i = -V_inf · (sin(alpha) - dz/dx_i · cos(alpha))
% (para perfil con curvatura; para NACA 0012 simétrico dz/dx ≠ 0 solo por espesor)

% -- Vectores planos para todos los paneles
Xv = zeros(N_tot, 1);
Xc = zeros(N_tot, 1);
Yc = zeros(N_tot, 1);
YL = zeros(N_tot, 1);
YR = zeros(N_tot, 1);
Zv = zeros(N_tot, 1);   % z del vórtice (extradós)
Zc = zeros(N_tot, 1);   % z del punto de control (extradós)
DZdx = zeros(N_tot, 1); % pendiente dz/dx en el punto de control
Dx = zeros(N_tot, 1);

for j = 1:N_y
    for i = 1:N_x
        k        = (j-1)*N_x + i;
        Xv(k)    = x_vort(i);
        Xc(k)    = x_ctrl(i);
        Yc(k)    = y_ctrl(j);
        YL(k)    = y_edges(j);
        YR(k)    = y_edges(j+1);
        Zv(k)    = z_vort(i);
        Zc(k)    = z_ctrl(i);
        DZdx(k)  = dzdx(i);
        Dx(k)    = dx_panel(i);
    end
end

%%
% =========================================================================
%  BLOQUE 3 — MATRIZ DE INFLUENCIA EN VUELO LIBRE (A_free)
%  Los puntos de control están en la superficie del extradós (Zc > 0)
%  La componente normal de la velocidad inducida incluye la inclinación local
% =========================================================================

fprintf('Calculando A_free (%d x %d)... ', N_tot, N_tot);
tic;
A_free = zeros(N_tot, N_tot);

for ks = 1:N_tot
    for kc = 1:N_tot
        % Punto de control en la superficie NACA 0012
        P = [Xc(kc), Yc(kc), Zc(kc)];
        V = horseshoe_velocity(P, Xv(ks), YL(ks), YR(ks), Zv(ks), L_wake);
        % Componente normal al perfil: n = (-dz/dx, 0, 1)/|n|
        % Proyección: w_n = (-dzdx * Vx + Vz) / norm(n)
        n_norm = sqrt(1 + DZdx(kc)^2);
        A_free(kc, ks) = (-DZdx(kc) * V(1) + V(3)) / n_norm;
    end
end
t_free = toc;
fprintf('OK (%.2f s)\n', t_free);

%%
% =========================================================================
%  BLOQUE 4 — TÉRMINO INDEPENDIENTE (RHS) CON PENDIENTE NACA 0012
%  Condición: (V_inf·n_freestream) + sum(Γ·A) = 0
%  V_inf_vector = V_inf*(cos α, 0, sin α)
%  n_panel = (-dz/dx, 0, 1)/|n|
%  RHS_k = -V_inf*(sin α - dz/dx_k * cos α) / |n_k|
% =========================================================================

RHS = zeros(N_tot, 1);
for k = 1:N_tot
    n_norm  = sqrt(1 + DZdx(k)^2);
    RHS(k)  = -V_inf * (sin(alpha) - DZdx(k) * cos(alpha)) / n_norm;
end

%%
% =========================================================================
%  BLOQUE 5 — VUELO LIBRE (referencia h → ∞)
% =========================================================================
Gamma_free  = A_free \ RHS;
[CL_free, CDi_free, CM_free, L_D_free] = compute_forces(Gamma_free, N_x, N_y, ...
                                               y_edges, Xv, rho, V_inf, ...
                                               alpha, S_ref, c, e_oswald);

fprintf('\n========== VUELO LIBRE (h → ∞) — NACA 0012 ==========\n');
fprintf('  CL      = %7.4f\n', CL_free);
fprintf('  CDi     = %7.5f  (%.2f counts)\n', CDi_free, CDi_free*1e4);
fprintf('  CM(c/4) = %+7.4f\n', CM_free);
fprintf('  L/D     = %7.2f\n', L_D_free);
fprintf('======================================================\n\n');

%%
% =========================================================================
%  BLOQUE 6 — LOOP PRINCIPAL DE ALTURAS
%  El suelo se sitúa en z = 0; el ala vuela a altura h sobre el suelo,
%  por lo que los puntos de control están en z = h + Zc(k).
%  La imagen especular del vórtice está en z = -(h + Zv(ks)).
% =========================================================================
CL_vec   = zeros(N_h, 1);
CDi_vec  = zeros(N_h, 1);
CM_vec   = zeros(N_h, 1);
L_N_vec  = zeros(N_h, 1);
L_D_vec  = zeros(N_h, 1);
Gamma_all = zeros(N_tot, N_h);

fprintf('  h/c   |   CL    | ΔCL%%   |  CDi×10³ |    CM    |   L/D   \n');
fprintf('--------|---------|--------|----------|----------|---------\n');

for ih = 1:N_h
    h = h_ratio(ih) * c;   % Altura del plano z=0 del perfil al suelo [m]

    A_img = zeros(N_tot, N_tot);
    for ks = 1:N_tot
        for kc = 1:N_tot
            % Punto de control a altura real sobre el suelo
            P_cp = [Xc(kc), Yc(kc), h + Zc(kc)];
            % Vórtice imagen: misma x,y pero z espejada con signo -
            z_img = -(h + Zv(ks));
            V_img = horseshoe_velocity(P_cp, Xv(ks), YL(ks), YR(ks), z_img, L_wake);
            % Proyección normal (misma normal que el panel real)
            n_norm = sqrt(1 + DZdx(kc)^2);
            A_img(kc, ks) = -(-DZdx(kc) * V_img(1) + V_img(3)) / n_norm;
        end
    end

    A_total = A_free + A_img;
    Gamma   = A_total \ RHS;
    Gamma_all(:, ih) = Gamma;

    [CL, CDi, CM, L_D] = compute_forces(Gamma, N_x, N_y, y_edges, Xv, ...
                                    rho, V_inf, alpha, S_ref, c, e_oswald);

    CL_vec(ih)  = CL;
    CDi_vec(ih) = CDi;
    CM_vec(ih)  = CM;
    L_D_vec(ih) = L_D;
    L_N_vec(ih) = CL * (0.5 * rho * V_inf^2 * S_ref);

    dCL_pct = (CL - CL_free) / CL_free * 100;
    fprintf(' %6.3f | %7.4f | %+5.1f%% | %8.4f | %+8.4f | %7.2f\n', ...
            h_ratio(ih), CL, dCL_pct, CDi*1e3, CM, L_D);
end
fprintf('\n');

%%
% =========================================================================
%  BLOQUE 7 — ANÁLISIS DE RESULTADOS
% =========================================================================
[dCL_max, idx_max] = max((CL_vec - CL_free) / CL_free * 100);
h_max = h_ratio(idx_max);

fprintf('========== RESUMEN EFECTO SUELO — NACA 0012 ==========\n');
fprintf('Incremento máximo de CL: %+.1f%% a h/c = %.2f\n', dCL_max, h_max);
fprintf('CL en h/c=%.2f: %.4f (vs %.4f libre)\n', h_ratio(1), CL_vec(1), CL_free);
fprintf('Reducción de CDi a h/c=%.2f: %.1f%%\n', ...
        h_ratio(1), (1 - CDi_vec(1)/CDi_free)*100);
fprintf('Mejora L/D a h/c=%.2f: %.1f%% (%.2f vs %.2f)\n', ...
        h_ratio(1), (L_D_vec(1)/L_D_free - 1)*100, L_D_vec(1), L_D_free);
fprintf('=======================================================\n\n');

%%
% =========================================================================
%  BLOQUE 8 — DISTRIBUCIÓN DE ΔCp EN LA ALTURA MÁS BAJA
% =========================================================================
Gamma_lowest = Gamma_all(:, 1);
DeltaCp      = zeros(N_x, N_y);
for j = 1:N_y
    for i = 1:N_x
        k = (j-1)*N_x + i;
        DeltaCp(i, j) = 2 * Gamma_lowest(k) / (V_inf * Dx(k));
    end
end

%%
% =========================================================================
%  BLOQUE 9 — GRÁFICAS
% =========================================================================

% ---- Figura 1: Coeficientes Aerodinámicos vs h/c ----
fig1 = figure('Color','w','Position',[50, 50, 1400, 900], ...
              'Name','VLM NACA 0012 — Efecto Suelo');

h_fine  = linspace(h_ratio(1), h_ratio(end), 400);
CL_fit  = pchip(h_ratio, CL_vec, h_fine);
CDi_fit = pchip(h_ratio, CDi_vec, h_fine);
CM_fit  = pchip(h_ratio, CM_vec, h_fine);
LD_fit  = pchip(h_ratio, L_D_vec, h_fine);

% 9.1 — CL vs h/c
ax1 = subplot(2,3,1);
fill([h_fine, fliplr(h_fine)], [CL_fit, CL_free*ones(size(h_fine))], ...
     [0.2, 0.6, 1.0], 'FaceAlpha', 0.18, 'EdgeColor','none'); hold on;
plot(h_fine, CL_fit, 'b-', 'LineWidth', 2.8);
plot(h_ratio, CL_vec, 'bo', 'MarkerSize', 9, 'MarkerFaceColor','b');
yline(CL_free, 'r--', 'LineWidth', 2.0, 'Label', sprintf('C_L^{\\infty} = %.3f', CL_free));
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('C_L', 'FontSize', 12, 'FontWeight','bold');
title('Coeficiente de Sustentación', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.2 — Incremento porcentual de CL
ax2 = subplot(2,3,2);
dCL_pct = (CL_vec - CL_free) / CL_free * 100;
dCL_fit = (CL_fit - CL_free) / CL_free * 100;
plot(h_fine, dCL_fit, 'Color', [0.8,0.3,0.3], 'LineWidth', 2.8); hold on;
plot(h_ratio, dCL_pct, 's', 'Color', [0.8,0.3,0.3], 'MarkerSize', 9, ...
     'MarkerFaceColor', [1,0.5,0.5]);
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('\DeltaC_L / C_L^{\infty} [%]', 'FontSize', 12, 'FontWeight','bold');
title('Incremento de Sustentación', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.3 — CDi vs h/c
ax3 = subplot(2,3,3);
plot(h_fine, CDi_fit*1e3, 'r-', 'LineWidth', 2.8); hold on;
%plot(h_ratio, CDi_vec*1e3, 'rs', 'MarkerSize', 9, 'MarkerFaceColor','r');
yline(CDi_free*1e3, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('C_{Di} \times 10^{3}', 'FontSize', 12, 'FontWeight','bold');
title('Drag Inducido', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.4 — CM vs h/c
ax4 = subplot(2,3,4);
plot(h_fine, CM_fit, 'Color', [0.2,0.7,0.3], 'LineWidth', 2.8); hold on;
%plot(h_ratio, CM_vec, '^', 'Color', [0.2,0.7,0.3], 'MarkerSize', 9,'MarkerFaceColor', [0.4,0.9,0.5]);
yline(CM_free, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('C_M (c/4)', 'FontSize', 12, 'FontWeight','bold');
title('Momento de Cabeceo', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.5 — L/D vs h/c
ax5 = subplot(2,3,5);
plot(h_fine, LD_fit, 'Color', [0.6,0.2,0.8], 'LineWidth', 2.8); hold on;
plot(h_ratio, L_D_vec, 'd', 'Color', [0.6,0.2,0.8], 'MarkerSize', 9, ...
     'MarkerFaceColor', [0.8,0.5,0.9]);
yline(L_D_free, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('L/D', 'FontSize', 12, 'FontWeight','bold');
title('Eficiencia Aerodinámica', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.6 — Distribución de circulación
ax6 = subplot(2,3,6);
colors_h = cool(N_h);
for ih = 1:min(6, N_h)
    Gamma_h = Gamma_all(:, ih);
    Gamma_strip_h = zeros(N_y, 1);
    for j = 1:N_y
        Gamma_strip_h(j) = sum(Gamma_h((j-1)*N_x+1 : j*N_x));
    end
    plot(y_ctrl/b, Gamma_strip_h, '-', 'Color', colors_h(ih,:), ...
         'LineWidth', 2.0, 'DisplayName', sprintf('h/c=%.2f', h_ratio(ih)));
    hold on;
end
Gf_strip = zeros(N_y, 1);
for j = 1:N_y
    Gf_strip(j) = sum(Gamma_free((j-1)*N_x+1 : j*N_x));
end
plot(y_ctrl/b, Gf_strip, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Libre');
xlabel('2y/b', 'FontSize', 12, 'FontWeight','bold');
ylabel('\Gamma_{strip} [m²/s]', 'FontSize', 12, 'FontWeight','bold');
title('Distribución de Circulación', 'FontSize', 13, 'FontWeight','bold');
legend('Location', 'north', 'NumColumns', 3, 'FontSize', 9);
grid on; box on; xlim([-1, 1]);

sgtitle(sprintf(['VLM NACA 0012 — Efecto Suelo   |   \\alpha = %.0f°   |   ' ...
                 'AR = %.1f   |   e = %.2f   |   N_y = %d'], ...
                rad2deg(alpha), AR, e_oswald, N_y), ...
        'FontSize', 15, 'FontWeight', 'bold');

% ---- Figura 2: Mapa 3D de ΔCp ----
fig2 = figure('Color','w','Position',[450, 100, 950, 650], ...
              'Name', sprintf('VLM NACA 0012 — Delta Cp (h/c=%.2f)', h_ratio(1)));

[Xmesh, Ymesh] = meshgrid(x_ctrl/c, y_ctrl);
surf(Xmesh', Ymesh', DeltaCp, 'EdgeColor','none', 'FaceColor','interp');
colormap(turbo);
cb = colorbar;
cb.Label.String = '\DeltaC_p';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

xlabel('x/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('y [m]', 'FontSize', 12, 'FontWeight','bold');
zlabel('\DeltaC_p', 'FontSize', 12, 'FontWeight','bold');
title(sprintf('Distribución de Presión   (h/c = %.2f,  \\alpha = %.0f°)', ...
              h_ratio(1), rad2deg(alpha)), 'FontSize', 14, 'FontWeight','bold');
view(-40, 30); grid on; axis tight;

hold on;
patch([-0.1, 1.1, 1.1, -0.1], [-b/2-0.3, -b/2-0.3, b/2+0.3, b/2+0.3], ...
      min(DeltaCp(:))*[1,1,1,1], 'FaceColor',[0.6,0.6,0.6], ...
      'FaceAlpha',0.3, 'EdgeColor','none');

% ---- Figura 3: Sustentación en Newtons ----
figure('Color','w','Position',[150, 150, 900, 550], 'Name', 'Fuerza de Sustentación — NACA 0012');

h_fine_N = linspace(h_ratio(1), h_ratio(end), 400);
L_N_fit  = pchip(h_ratio, L_N_vec, h_fine_N);
L_free_N = CL_free * (0.5 * rho * V_inf^2 * S_ref);

fill([h_fine_N, fliplr(h_fine_N)], [L_N_fit, L_free_N*ones(size(h_fine_N))], ...
     [1.0, 0.5, 0.3], 'FaceAlpha', 0.15, 'EdgeColor','none'); hold on;

plot(h_fine_N, L_N_fit, '-', 'Color', [0.8,0.2,0.1], 'LineWidth', 3.5);
yline(L_free_N, 'k--', 'LineWidth', 2.2, 'Label', 'Vuelo Libre', ...
      'LabelHorizontalAlignment','left');

title('Fuerza de Sustentación vs Altura — NACA 0012', 'FontSize', 15, 'FontWeight','bold');
xlabel('Altura relativa (h/c)', 'FontSize', 13, 'FontWeight','bold');
ylabel('Sustentación L [N]', 'FontSize', 13, 'FontWeight','bold');
grid on; box on;



% % ---- Figura 4: Geometría 3D del ala NACA 0012 ----
% fig4 = figure('Color','w','Position',[200, 200, 1000, 600], ...
%               'Name', 'Geometría 3D — Ala NACA 0012');
% hold on; grid on; view(3);
% 
% % Dibujar paneles del extradós e intradós
% N_x_vis = 10; N_y_vis = 20;
% beta_cv  = linspace(0, pi, N_x_vis + 1);
% x_ev     = c/2 * (1 - cos(beta_cv));
% beta_sv  = linspace(0, pi, N_y_vis + 1);
% y_ev     = -b/2 * cos(beta_sv);
% 
% for j = 1:N_y_vis
%     for i = 1:N_x_vis
%         x1=x_ev(i); x2=x_ev(i+1);
%         y1=y_ev(j); y2=y_ev(j+1);
%         z11=naca_z(x1,c,0.12); z21=naca_z(x2,c,0.12);
%         % Extradós
%         fill3([x1 x2 x2 x1],[y1 y1 y2 y2],[z11 z21 z21 z11], ...
%               [0.75 0.85 0.95],'EdgeColor','k','FaceAlpha',0.85);
%         % Intradós
%         fill3([x1 x2 x2 x1],[y1 y1 y2 y2],[-z11 -z21 -z21 -z11], ...
%               [0.65 0.75 0.85],'EdgeColor','k','FaceAlpha',0.70);
%     end
% end
% 
% % Perfiles en sección central y tips
% xi_p = linspace(0, c, 300);
% zp   = arrayfun(@(x) naca_z(x,c,0.12), xi_p);
% for yp = [0, -b/2, b/2]
%     col = [0.1 0.3 0.7]; if abs(yp)>0, col=[0.7 0.1 0.1]; end
%     plot3(xi_p, yp*ones(size(xi_p)),  zp, 'Color', col, 'LineWidth', 1.8);
%     plot3(xi_p, yp*ones(size(xi_p)), -zp, 'Color', col, 'LineWidth', 1.8);
% end
% 
% xlabel('x [m]','FontSize',11,'FontWeight','bold');
% ylabel('y [m]','FontSize',11,'FontWeight','bold');
% zlabel('z [m]','FontSize',11,'FontWeight','bold');
% title('Geometría 3D — Ala NACA 0012 (VLM)','FontSize',13,'FontWeight','bold');
% axis equal; xlim([0 c]); ylim([-b/2 b/2]); zlim([-0.15*c 0.15*c]);

%%
% =========================================================================
%  BLOQUE 10 — EXPORTAR TABLA FINAL
% =========================================================================
fprintf('\n========== TABLA COMPLETA DE RESULTADOS — NACA 0012 ==========\n');
fprintf('%6s | %8s | %10s | %10s | %8s | %8s\n', ...
        'h/c', 'CL', 'CDi', 'CM', 'L/D', 'ΔCL[%]');
fprintf('%s\n', repmat('-',1,68));
for ih = 1:N_h
    dCL = (CL_vec(ih) - CL_free) / CL_free * 100;
    fprintf('%6.3f | %8.4f | %10.6f | %+10.4f | %8.2f | %+7.1f\n', ...
            h_ratio(ih), CL_vec(ih), CDi_vec(ih), CM_vec(ih), L_D_vec(ih), dCL);
end
fprintf('%6s | %8.4f | %10.6f | %+10.4f | %8.2f | %7s\n', ...
        '∞', CL_free, CDi_free, CM_free, L_D_free, '0.0');
fprintf('==============================================================\n');

%%
% =========================================================================
%  FUNCIONES LOCALES
% =========================================================================
function [AR , e_oswald] = Oswald(b,c)
    
    AR       = b / c;      % Alargamiento

    % Factor de Oswald para AR bajo (empírico, validado experimentalmente)
    % Fórmula de Raymer: e ≈ 1.78*(1 - 0.045*AR^0.68) - 0.64
    if AR < 4
        e_oswald = 1.78 * (1 - 0.045 * AR^0.68) - 0.64;
    else
        e_oswald = 0.95;  % AR alto
    end
end

function z = naca_z(x_abs, c, t_ratio)
% Ordenada z del extradós del perfil NACA 4 dígitos simétrico
% Coeficientes estándar NACA, trailing edge cerrado (a4 = -0.1015)
    xi = max(0, min(1, x_abs / c));
    a0 =  0.2969;  a1 = -0.1260;
    a2 = -0.3516;  a3 =  0.2843;  a4 = -0.1015;
    z  = 5 * t_ratio * c * (a0*sqrt(xi) + a1*xi + a2*xi^2 + a3*xi^3 + a4*xi^4);
end

function dzx = naca_dzdx(x_abs, c, t_ratio)
% Pendiente dz/dx del extradós NACA 4 dígitos simétrico
    xi  = max(1e-6, min(1 - 1e-6, x_abs / c));
    a0 =  0.2969;  a1 = -0.1260;
    a2 = -0.3516;  a3 =  0.2843;  a4 = -0.1015;
    % dz/dxi * (1/c) = dz/dx
    dzx = (5 * t_ratio / c) * ...
          (a0 / (2*sqrt(xi)) + a1 + 2*a2*xi + 3*a3*xi^2 + 4*a4*xi^3);
end

function V = horseshoe_velocity(P, xv, yL, yR, zw, L_wake)
% Velocidad inducida en P por vórtice de herradura unitario (Γ=1)
    A     = [xv,          yL, zw];
    B     = [xv,          yR, zw];
    A_far = [xv + L_wake, yL, zw];
    B_far = [xv + L_wake, yR, zw];
    V  = biot_savart(P, A, B);
    V  = V + biot_savart(P, B, B_far);
    V  = V + biot_savart(P, A_far, A);
end

function V = biot_savart(P, A, B)
% Ley de Biot-Savart para segmento A→B con Γ=1
    EPSILON = 1e-10;
    r1  = P - A;   r2  = P - B;   r0  = B - A;
    cr  = cross(r1, r2);
    cr2 = dot(cr, cr);
    if cr2 < EPSILON^2, V = [0,0,0]; return; end
    r1n = norm(r1);   r2n = norm(r2);
    if r1n < EPSILON || r2n < EPSILON, V = [0,0,0]; return; end
    K = (1/(4*pi)) * (dot(r0,r1)/r1n - dot(r0,r2)/r2n) / cr2;
    V = K * cr;
end

function [CL, CDi, CM, L_D] = compute_forces(Gamma, N_x, N_y, y_edges, Xv, ...
                                              rho, V_inf, alpha, S_ref, c, e_oswald)
    L_total = 0;   M_total = 0;
    x_ref   = 0.25 * c;
    for j = 1:N_y
        dy = y_edges(j+1) - y_edges(j);
        for i = 1:N_x
            k  = (j-1)*N_x + i;
            dL = rho * V_inf * Gamma(k) * dy;
            L_total = L_total + dL;
            M_total = M_total + dL * (x_ref - Xv(k));
        end
    end
    b_span = y_edges(end) - y_edges(1);
    AR_loc = b_span^2 / S_ref;
    CL  = L_total / (0.5 * rho * V_inf^2 * S_ref);
    CDi = CL^2 / (pi * e_oswald * AR_loc);
    CM  = M_total / (0.5 * rho * V_inf^2 * S_ref * c);
    if CDi > 1e-8
        L_D = CL / CDi;
    else
        L_D = inf;
    end
end
