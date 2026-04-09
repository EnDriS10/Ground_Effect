% =========================================================================
%  Vortex Lattice Method (VLM) con Efecto Suelo 
% =========================================================================

clear; clc; close all;

%%
% =========================================================================
%  BLOQUE 1 — PARÁMETROS Y MALLA
% =========================================================================

% Geometría
b        = 10.0;        % Envergadura [m]
c        = 1.0;        % Cuerda [m]

% Vuelo
V_inf    = 15;         % Velocidad de la corriente libre [m/s]
alpha    = deg2rad(5); % Ángulo de ataque [rad]
rho      = 1.225;      % Densidad del aire [kg/m³]
S_ref    = b * c;      % Área de referencia [m²]

% Discretización (mayor resolución para AR bajo)
N_x      = 10;         % Paneles en dirección de la cuerda
N_y      = 48;         % Paneles en dirección de la envergadura (aumentado)
N_tot    = N_x * N_y;

% Estela
L_wake   = 150 * c; 

% Alturas sobre el suelo h/c
h_ratio  = [0.10, 0.15, 0.20, 0.25, 0.30,0.35,0.40,0.45, 0.50,0.55, ... 
    0.60,0.65,0.70,0.75,0.80,0.85,0.90, 1.00,1.25, 1.50,1.75, 2.00,2.50,3.00];

N_h = length(h_ratio);

[AR , e_oswald] = Oswald(b,c);
fprintf('Factor de Oswald calculado: e = %.3f (AR=%.1f)\n', e_oswald, AR);

fprintf('VLM Efecto Suelo | AR=%.1f | Nx=%d Ny=%d | Ntot=%d\n', ...
        AR, N_x, N_y, N_tot);

fprintf('Estela: %.0f×c | Oswald: e=%.3f\n\n', L_wake/c, e_oswald);

%Inicializacion de Malla don Distribucion Coseno
[Xv, Xc, Yc, YL, YR, Dx, x_edges, y_edges, x_ctrl, y_ctrl, dx_panel] = Malla_VLM(b, c, N_x, N_y);


%Término Independiente
RHS = -V_inf * sin(alpha) * ones(N_tot, 1);

%%
% =========================================================================
%  BLOQUE 2 — VUELO LIBRE
% =========================================================================

[A_free, t_free] = Matriz_free(N_tot, Xc, Yc, Xv, YL, YR, L_wake); %Matriz de influencia

Gamma_free  = A_free \ RHS;

[CL_free, CDi_free, CM_free, L_D_free] = compute_forces(Gamma_free, N_x, N_y, ...
                                               y_edges, Xv, rho, V_inf, ...
                                               alpha, S_ref, c, e_oswald);

fprintf('\n========== VUELO LIBRE (h → ∞) ==========\n');
fprintf('  CL      = %7.4f\n', CL_free);
fprintf('  CDi     = %7.5f  (%.2f counts)\n', CDi_free, CDi_free*1e4);
fprintf('  CM(c/4) = %+7.4f\n', CM_free);
fprintf('  L/D     = %7.2f\n', L_D_free);
fprintf('=========================================\n\n');


%%
% =========================================================================
%  BLOQUE 3 — LOOP PRINCIPAL DE ALTURAS
% =========================================================================

[CL_vec, CDi_vec, CM_vec, L_N_vec, L_D_vec, Gamma_all] = Main_loop(...
    h_ratio, c, N_h, N_tot, Xc, Yc, Xv, YL, YR, L_wake, A_free, RHS, ...
    N_x, N_y, y_edges, rho, V_inf, alpha, S_ref, e_oswald, CL_free);

% DISTRIBUCIÓN DE ΔCp
N_punto_altura=1;
Gamma_especifica = Gamma_all(:, N_punto_altura);  % h/c elegido
DeltaCp = compute_delta_cp(Gamma_especifica, N_x, N_y, V_inf, Dx);


% Análisis
% Incrementos máximos
[dCL_max, idx_max] = max((CL_vec - CL_free) / CL_free * 100);
h_max = h_ratio(idx_max);

fprintf('========== RESUMEN EFECTO SUELO ==========\n');
fprintf('Incremento máximo de CL: %+.1f%% a h/c = %.2f\n', dCL_max, h_max);
fprintf('CL en h/c=%.2f: %.4f (vs %.4f libre)\n', h_ratio(1), CL_vec(1), CL_free);
fprintf('Reducción de CDi a h/c=%.2f: %.1f%%\n', ...
        h_ratio(1), (1 - CDi_vec(1)/CDi_free)*100);
fprintf('Mejora L/D a h/c=%.2f: %.1f%% (%.2f vs %.2f)\n', ...
        h_ratio(1), (L_D_vec(1)/L_D_free - 1)*100, L_D_vec(1), L_D_free);
fprintf('==========================================\n\n');

%%
% =========================================================================
%  BLOQUE 4 — GRÁFICAS
% =========================================================================


% ---- Figura 1: Coeficientes Aerodinámicos vs h/c ----
fig1 = figure('Color','w','Position',[50, 50, 1400, 900], ...
              'Name','VLM — Efecto Suelo (Mejorado)');

% Ajuste suave para curvas
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
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('\DeltaC_L / C_L^{\infty} [%]', 'FontSize', 12, 'FontWeight','bold');
title('Incremento de Sustentación', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.3 — CDi vs h/c
ax3 = subplot(2,3,3);
plot(h_fine, CDi_fit*1e3, 'r-', 'LineWidth', 2.8); hold on;
yline(CDi_free*1e3, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('C_{Di} \times 10^{3}', 'FontSize', 12, 'FontWeight','bold');
title('Drag Inducido', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.4 — CM vs h/c
ax4 = subplot(2,3,4);
plot(h_fine, CM_fit, 'Color', [0.2,0.7,0.3], 'LineWidth', 2.8); hold on;
yline(CM_free, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('C_M (c/4)', 'FontSize', 12, 'FontWeight','bold');
title('Momento de Cabeceo', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.5 — L/D vs h/c
ax5 = subplot(2,3,5);
plot(h_fine, LD_fit, 'Color', [0.6,0.2,0.8], 'LineWidth', 2.8); hold on;
yline(L_D_free, 'k--', 'LineWidth', 1.8, 'Label', 'Vuelo libre');
xlabel('h/c', 'FontSize', 12, 'FontWeight','bold');
ylabel('L/D', 'FontSize', 12, 'FontWeight','bold');
title('Eficiencia Aerodinámica', 'FontSize', 13, 'FontWeight','bold');
grid on; box on; xlim([0, h_ratio(end)+0.2]);

% 9.6 — Distribución de circulación
ax6 = subplot(2,3,6);
colors_h = cool(N_h);
for ih = 1:min(6, N_h)  % Solo primeras 6 alturas para claridad
    Gamma_h = Gamma_all(:, ih);
    Gamma_strip_h = zeros(N_y, 1);
    for j = 1:N_y
        Gamma_strip_h(j) = sum(Gamma_h((j-1)*N_x+1 : j*N_x));
    end
    plot(y_ctrl/b, Gamma_strip_h, '-', 'Color', colors_h(ih,:), ...
         'LineWidth', 2.0, 'DisplayName', sprintf('h/c=%.2f', h_ratio(ih)));
    hold on;
end
% Vuelo libre
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

sgtitle(sprintf(['VLM Mejorado — Efecto Suelo   |   \\alpha = %.0f°   |   ' ...
                 'AR = %.1f   |   e = %.2f   |   N_y = %d'], ...
                rad2deg(alpha), AR, e_oswald, N_y), ...
        'FontSize', 15, 'FontWeight', 'bold');

% ---- Figura 2: Mapa 3D de ΔCp ----
fig2 = figure('Color','w','Position',[450, 100, 950, 650], ...
              'Name', sprintf('VLM — Delta Cp (h/c=%.2f)', h_ratio(N_punto_altura)));

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

% Representación del suelo
hold on;
patch([-0.1, 1.1, 1.1, -0.1], [-b/2-0.3, -b/2-0.3, b/2+0.3, b/2+0.3], ...
      min(DeltaCp(:))*[1,1,1,1], 'FaceColor',[0.6,0.6,0.6], ...
      'FaceAlpha',0.3, 'EdgeColor','none');

% ---- Figura 3: Sustentación en Newtons ----
figure('Color','w','Position',[150, 150, 900, 550], 'Name', 'Fuerza de Sustentación');

h_fine_N = linspace(h_ratio(1), h_ratio(end), 400);
L_N_fit = pchip(h_ratio, L_N_vec, h_fine_N);
L_free_N = CL_free * (0.5 * rho * V_inf^2 * S_ref);

fill([h_fine_N, fliplr(h_fine_N)], [L_N_fit, L_free_N*ones(size(h_fine_N))], ...
     [1.0, 0.5, 0.3], 'FaceAlpha', 0.15, 'EdgeColor','none'); hold on;

plot(h_fine_N, L_N_fit, '-', 'Color', [0.8,0.2,0.1], 'LineWidth', 3.5);
yline(L_free_N, 'k--', 'LineWidth', 2.2, 'Label', 'Vuelo Libre', ...
      'LabelHorizontalAlignment','left');

title('Fuerza de Sustentación vs Altura', 'FontSize', 15, 'FontWeight','bold');
xlabel('Altura relativa (h/c)', 'FontSize', 13, 'FontWeight','bold');
ylabel('Sustentación L [N]', 'FontSize', 13, 'FontWeight','bold');
grid on; box on;



%%
% =========================================================================
%  BLOQUE 5 — EXPORTAR TABLA FINAL
% =========================================================================
fprintf('\n========== TABLA COMPLETA DE RESULTADOS ==========\n');
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
fprintf('==================================================\n');

%%

% =========================================================================
%  BLOQUE 6 — FUNCIONES LOCALES
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


function [Xv, Xc, Yc, YL, YR, Dx, x_edges, y_edges, x_ctrl, y_ctrl, dx_panel] = Malla_VLM(b, c, N_x, N_y)

% GENERAR_MALLA_VLM: Discretiza el ala en paneles usando una distribución coseno.

    N_tot = N_x * N_y;

    % -- Cuerda: distribución coseno (mayor densidad en LE y TE)
    beta_c   = linspace(0, pi, N_x + 1);
    x_edges  = c/2 * (1 - cos(beta_c));          % [0 → c]
    dx_panel = diff(x_edges);                    % Longitud de cada panel en x
    x_vort   = x_edges(1:N_x) + dx_panel / 4;    % 1/4 del panel (vórtice)
    x_ctrl   = x_edges(1:N_x) + 3*dx_panel / 4;  % 3/4 del panel (control)

    % -- Envergadura: distribución coseno (mejor resolución en los wingtips)
    beta_s   = linspace(0, pi, N_y + 1);
    y_edges  = -b/2 * cos(beta_s);               % [-b/2 → +b/2]
    y_ctrl   = 0.5 * (y_edges(1:N_y) + y_edges(2:N_y+1));

    % -- Vectores para todos los paneles (almacenamiento plano)
    Xv = zeros(N_tot, 1);   % x del vórtice ligado
    Xc = zeros(N_tot, 1);   % x del punto de control
    Yc = zeros(N_tot, 1);   % y del punto de control
    YL = zeros(N_tot, 1);   % y borde izquierdo del panel (envergadura)
    YR = zeros(N_tot, 1);   % y borde derecho del panel (envergadura)
    Dx = zeros(N_tot, 1);   % Longitud del panel en cuerda

    for j = 1:N_y
        for i = 1:N_x
            k        = (j-1)*N_x + i;
            Xv(k)    = x_vort(i);
            Xc(k)    = x_ctrl(i);
            Yc(k)    = y_ctrl(j);
            YL(k)    = y_edges(j);
            YR(k)    = y_edges(j+1);
            Dx(k)    = dx_panel(i);
        end
    end
end

function [A_free, t_free] = Matriz_free(N_tot, Xc, Yc, Xv, YL, YR, L_wake)

%   A_free - Matriz de coeficientes de influencia aerodinámica [N_tot x N_tot]
%   t_free - Tiempo de ejecución en segundos

    fprintf('Calculando A_free (%d x %d)... ', N_tot, N_tot);
    tic;
    
    A_free = zeros(N_tot, N_tot);
    
    for ks = 1:N_tot          % Panel fuente (herradura)
        for kc = 1:N_tot      % Panel de control (punto de evaluación)
            P = [Xc(kc), Yc(kc), 0.0];
            
            % Calcula la velocidad inducida usando la función auxiliar
            V = horseshoe_velocity(P, Xv(ks), YL(ks), YR(ks), 0.0, L_wake);
            
            % Almacena solo la componente normal (z)
            A_free(kc, ks) = V(3);   
        end
    end
    
    t_free = toc;
    fprintf('OK (%.2f s)\n', t_free);
end

function V = horseshoe_velocity(P, xv, yL, yR, zw, L_wake)
% -------------------------------------------------------------------------
%  Velocidad inducida en P por vórtice de herradura unitario (Γ=1)
% -------------------------------------------------------------------------
    A     = [xv,          yL, zw];
    B     = [xv,          yR, zw];
    A_far = [xv + L_wake, yL, zw];
    B_far = [xv + L_wake, yR, zw];

    V  = compute_circulation(P, A, B);        % Vórtice ligado: A → B
    V  = V + compute_circulation(P, B, B_far); % Pierna derecha: B → B_far
    V  = V + compute_circulation(P, A_far, A); % Pierna izquierda: A_far → A
end


function V = compute_circulation(P, A, B)
% -------------------------------------------------------------------------
%  Ley de Biot-Savart para segmento A→B con Γ=1
%  Regularización mejorada para evitar singularidades
% -------------------------------------------------------------------------
    EPSILON = 1e-10;   % Núcleo del vórtice (regularización)

    r1  = P - A;
    r2  = P - B;
    r0  = B - A;

    cr  = cross(r1, r2);
    cr2 = dot(cr, cr);

    if cr2 < EPSILON^2
        V = [0, 0, 0];
        return;
    end

    r1n = norm(r1);
    r2n = norm(r2);

    if r1n < EPSILON || r2n < EPSILON
        V = [0, 0, 0];
        return;
    end

    K = (1 / (4*pi)) * (dot(r0, r1)/r1n - dot(r0, r2)/r2n) / cr2;
    V = K * cr;
end


function [CL, CDi, CM, L_D] = compute_forces(Gamma, N_x, N_y, y_edges, Xv, ...
                                              rho, V_inf, alpha, S_ref, c, e_oswald)

    L_total = 0;
    M_total = 0;
    x_ref   = 0.25 * c;   % Eje de momentos (c/4)

    for j = 1:N_y
        dy = y_edges(j+1) - y_edges(j);
        for i = 1:N_x
            k  = (j-1)*N_x + i;
            dL = rho * V_inf * Gamma(k) * dy;   % Kutta-Joukowski
            L_total = L_total + dL;
            M_total = M_total + dL * (x_ref - Xv(k));
        end
    end

    b   = y_edges(end) - y_edges(1);
    AR  = b^2 / S_ref;

    CL  = L_total / (0.5 * rho * V_inf^2 * S_ref);
    
    % CORRECCIÓN: Factor de Oswald e ≠ 1.0 (crucial para AR bajo)
    CDi = CL^2 / (pi * e_oswald * AR);
    
    CM  = M_total / (0.5 * rho * V_inf^2 * S_ref * c);
    
    % Eficiencia L/D (solo drag inducido, sin drag parasitic)
    if CDi > 1e-8
        L_D = CL / CDi;
    else
        L_D = inf;
    end
end

function [CL_vec, CDi_vec, CM_vec, L_N_vec, L_D_vec, Gamma_all] = Main_loop(...
    h_ratio, c, N_h, N_tot, Xc, Yc, Xv, YL, YR, L_wake, A_free, RHS, ...
    N_x, N_y, y_edges, rho, V_inf, alpha, S_ref, e_oswald, CL_free)

% EJECUTAR_LOOP_ALTURAS: Calcula el efecto suelo para diferentes alturas (h/c)
% mediante el método de las imágenes.

    % Inicialización de vectores de resultados
    CL_vec    = zeros(N_h, 1);
    CDi_vec   = zeros(N_h, 1);
    CM_vec    = zeros(N_h, 1);
    L_N_vec   = zeros(N_h, 1);
    L_D_vec   = zeros(N_h, 1);
    Gamma_all = zeros(N_tot, N_h);

    fprintf('  h/c   |   CL    | ΔCL%%   |  CDi×10³ |    CM    |   L/D   \n');
    fprintf('--------|---------|--------|----------|----------|---------\n');

    for ih = 1:N_h
        h = h_ratio(ih) * c;   % Altura física en metros
        
        % -- Matriz de influencia de la imagen --
        % Se evalúa la influencia de los vórtices imagen (en z = -h) 
        % sobre los puntos de control del ala real (en z = h).
        A_img = zeros(N_tot, N_tot);
        for ks = 1:N_tot
            for kc = 1:N_tot
                P_cp = [Xc(kc), Yc(kc), h];
                % El vórtice imagen está en -h
                V_img = horseshoe_velocity(P_cp, Xv(ks), YL(ks), YR(ks), -h, L_wake);
                % Condición de contorno: Gamma_img = -Gamma_real
                A_img(kc, ks) = -V_img(3); 
            end
        end

        % Resolución del sistema: (A_real + A_imagen) * Gamma = RHS
        A_total = A_free + A_img;
        Gamma   = A_total \ RHS;
        Gamma_all(:, ih) = Gamma;

        % Cálculo de fuerzas para la altura actual
        [CL, CDi, CM, L_D] = compute_forces(Gamma, N_x, N_y, y_edges, Xv, ...
                                        rho, V_inf, alpha, S_ref, c, e_oswald);
        
        % Almacenamiento en vectores
        CL_vec(ih)  = CL;
        CDi_vec(ih) = CDi;
        CM_vec(ih)  = CM;
        L_D_vec(ih) = L_D;
        
        % Sustentación en Newtons
        L_N_vec(ih) = CL * (0.5 * rho * V_inf^2 * S_ref); 
        
        % Cálculo de incremento para impresión en consola
        dCL_pct = (CL - CL_free) / CL_free * 100;
        fprintf(' %6.3f | %7.4f | %+5.1f%% | %8.4f | %+8.4f | %7.2f\n', ...
                h_ratio(ih), CL, dCL_pct, CDi*1e3, CM, L_D);
    end
    fprintf('\n');
end

function DeltaCp = compute_delta_cp(Gamma_especifica, N_x, N_y, V_inf, Dx)

    DeltaCp = zeros(N_x, N_y);

    for j = 1:N_y
        for i = 1:N_x
            % Índice lineal
            k = (j-1)*N_x + i;
            
            % Cálculo del coeficiente de presión diferencial
            DeltaCp(i, j) = 2 * Gamma_especifica(k) / (V_inf * Dx(k));
        end
    end
end