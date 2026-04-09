

clear; clc; close all;
figure('Color', 'k', 'Position', [100 100 1000 700]);

% 1. Definir variables
b = 10;           % Envergadura (m)
Gamma = 25;       % Circulación (m^2/s)
V_viento = 15;       % Velocidad del viento (m/s)
L_estela = 80;    % Largo estela (m)
rho = 1.225;      % Densidad aire (kg/m^3)
h = 4;            % Altura del ala sobre el plano (m)

% 2. Solución analítica
w_inf = Gamma / (-pi * b); 
factor_h = (16 * (h/b)^2) / (1 + 16 * (h/b)^2);
w_total_analitico = w_inf * factor_h;

L_inf = rho * V_viento * Gamma * b;
L_analitico = L_inf * (w_inf / w_total_analitico);

fprintf('----- Resultados Analíticos: -----\n');
fprintf('Sustentación en Aire Libre (L_inf):  %.2f N\n', L_inf);
fprintf('Downwash Analítico Efecto Suelo:       %.4f m/s\n', w_total_analitico);
fprintf('Sustentación Analítica (L_analitico):      %.2f N\n', L_analitico);

% 3. Creamos malla tridimensional
[X, Y, Z] = meshgrid(linspace(-10, 40, 35), ...
                     linspace(-b*1.2, b*1.2, 30), ...
                     linspace(-h*2.2, h*2.2, 30));

U = V_viento * ones(size(X)); V = zeros(size(X)); W = zeros(size(X));        

% 4. Aplicamos m. imágenes
real_fil = {[L_estela, -b/2, h; 0, -b/2, h], [0, -b/2, h; 0, b/2, h], [0, b/2, h; L_estela, b/2, h]};
espj_fil = {[L_estela, -b/2, -h; 0, -b/2, -h], [0, -b/2, -h; 0, b/2, -h], [0, b/2, -h; L_estela, b/2, -h]};

filamentos = [real_fil, espj_fil];
gammas = [Gamma, Gamma, Gamma, -Gamma, -Gamma, -Gamma];

fprintf('\n Aplicando Biot-Savart... espere.\n');

for f = 1:6
    seg = filamentos{f}; G = gammas(f);
    A = seg(1,:); B = seg(2,:);
    for i = 1:numel(X)
        P = [X(i), Y(i), Z(i)];
        r1 = P - A; r2 = P - B; r0 = B - A;
        c = cross(r1, r2); mag_c2 = sum(c.^2);
        if mag_c2 > 0.1 
            V_ind = (G / (4 * pi)) * (c / mag_c2) * (dot(r1, r0)/norm(r1) - dot(r2, r0)/norm(r2));
            U(i) = U(i) + V_ind(1); V(i) = V(i) + V_ind(2); W(i) = W(i) + V_ind(3);
        end
    end
end

% 5. Sacamos L_numerico através del campo de velocidades en x=30

L_numerico = calcular_L_numerico(h, b, Gamma, V_viento, rho, L_estela);

w_numerico = w_inf / (L_numerico / L_inf);

error_numerico = 100 * abs((L_analitico - L_numerico) / L_analitico);

fprintf('\n----- Resultados Numéricos -----\n');
fprintf('Downwash Simulado Exacto:            %.4f m/s\n', w_numerico);
fprintf('Sustentación Simulada Exacta:        %.2f N\n', L_numerico);
fprintf('Error REAL Teórico vs Simulado:      %.4f%%\n\n', error_numerico);

%  6. Visualizacion 3d
Vel_mag = sqrt(U.^2 + V.^2 + W.^2);
hold on;

z_seeds = linspace(-h*1.8, h*1.8, 14); 
y_seeds = linspace(-b, b, 12);
[seedY, seedZ] = meshgrid(y_seeds, z_seeds);
seedX = -10 * ones(size(seedY));

xyz = stream3(X, Y, Z, U, V, W, seedX, seedY, seedZ);
colormap(turbo);

for i = 1:length(xyz)
    pts = xyz{i}; if isempty(pts), continue; end
    c_data = interp3(X, Y, Z, Vel_mag, pts(:,1), pts(:,2), pts(:,3));
    surface([pts(:,1), pts(:,1)], [pts(:,2), pts(:,2)], [pts(:,3), pts(:,3)], ...
            [c_data, c_data], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 0.8);
end

% Dibujando el ala fina
fill3([0 0 0.5 0.5], [-b/2 b/2 b/2 -b/2], [h h h h], [0.6 0.6 0.6], 'EdgeColor', 'w'); % Real
fill3([0 0 0.5 0.5], [-b/2 b/2 b/2 -b/2], [-h -h -h -h], [0.4 0.1 0.1], 'EdgeColor', 'r'); % Espejo

%Colores y encuadres
axis equal tight; view([-55, 20]); grid on;
cb = colorbar('Color', 'w', 'Location', 'eastoutside');
ylabel(cb, 'Velocidad Total (V_{inf} + V_{ind}) [m/s]', 'FontSize', 11);
clim([V_viento*0.9, V_viento*1.5]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
title(['Efecto Suelo: Lift Teórico = ', num2str(round(L_analitico)), ' N | Lift Simulado = ', num2str(round(L_numerico)), ' N'], 'Color', 'w');
%}

% =========================================================================
% ----- Extensión dinámica (Heun) (F_resultante= L-mg+(F_roz) -----
% =========================================================================
fprintf('\n--- SIMULANDO DINÁMICA DE VUELO (RK2 Numérico) ---\n');
% Variables extra
m = 800;          
g = 9.81;         
p = 600;          
dt = 0.05;        
t_final = 25;     
t = 0:dt:t_final; 

Y = zeros(2, length(t)); 
Y(:, 1) = [15; 0];  

for i = 1:(length(t)-1)
    h_i = Y(1, i);
    v_i = Y(2, i);
    
    % --- K1: Instante actual usando el SOLUCIONADOR NUMÉRICO ---
    % Descartar h=0 por el bien de la trama
    h_eval = max(h_i, 0.1); 
    L1 = calcular_L_numerico(h_eval, b, Gamma, V_viento, rho, L_estela);    
    a1 = (L1 - m*g - p*v_i) / m; 
    
    k1_h = v_i;                
    k1_v = a1;                 
    
    % --- Paso intermedio (Punto Medio) ---
    h_mid = h_i + k1_h * (dt / 2);
    v_mid = v_i + k1_v * (dt / 2);
    
    % --- K2: Punto medio usando el SOLUCIONADOR NUMÉRICO ---
    h_mid_eval = max(h_mid, 0.1);
    L2 = calcular_L_numerico(h_mid_eval, b, Gamma, V_viento, rho, L_estela);    
    a2 = (L2 - m*g - p*v_mid) / m;
    
    k2_h = v_mid;
    k2_v = a2;
    
    
    Y(1, i+1) = h_i + k2_h * dt;
    Y(2, i+1) = v_i + k2_v * dt;
    

    if Y(1, i+1) < 0.2
        Y(1, i+1) = 0.2;
        Y(2, i+1) = -0.1 * Y(2, i+1); 
    end
end

%Gráfica
figure('Color', 'w', 'Position', [200, 200, 800, 400]);
hold on;


p1 = plot(t, Y(1, :), 'b-', 'LineWidth', 2.5); 


factor_eq = L_inf / (m * g);
if factor_eq < 1 
    h_eq = sqrt(factor_eq / (16 - 16*factor_eq)) * b;
   
    p2 = yline(h_eq, 'r--', 'LineWidth', 2);
end


title(['Dinámica de Vuelo: Amortiguamiento en Efecto Suelo (p = ', num2str(p), ')'], 'FontSize', 14,'Color','k');
xlabel('Tiempo, t (s)', 'FontSize', 12, 'Color', 'k');
ylabel('Altura sobre el suelo, h (m)', 'FontSize', 12, 'Color', 'k');


set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
         'FontSize', 11, 'GridAlpha', 0.15, 'LineWidth', 1);
grid on;
box on; % 


xlim([min(t), max(t)]);

%Poniendo bonita la leyenda
if factor_eq < 1
    lgd = legend([p1, p2], 'Trayectoria Simulada (RK2)', 'Altura de Equilibrio Teórica', 'Location', 'northeast', 'FontSize', 11);
else
    lgd = legend(p1, 'Trayectoria Simulada (RK2)', 'Location', 'northeast', 'FontSize', 11);
end
lgd.Color = 'w';
lgd.TextColor = 'k';
lgd.EdgeColor = 'k';

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L_num = calcular_L_numerico(h_actual, b, Gamma, V_inf, rho, L_estela)
    
    w_inf = Gamma / (-pi * b); 
    L_inf = rho * V_inf * Gamma * b;
    
    % M. Imagen
    real_fil = {[L_estela, -b/2, h_actual; 0, -b/2, h_actual], ...
                [0, -b/2, h_actual; 0, b/2, h_actual], ...
                [0, b/2, h_actual; L_estela, b/2, h_actual]};
    espj_fil = {[L_estela, -b/2, -h_actual; 0, -b/2, -h_actual], ...
                [0, -b/2, -h_actual; 0, b/2, -h_actual], ...
                [0, b/2, -h_actual; L_estela, b/2, -h_actual]};
                
    filamentos = [real_fil, espj_fil];
    gammas = [Gamma, Gamma, Gamma, -Gamma, -Gamma, -Gamma];
    
    % 3. Punto de medida (X=40, Y=0, Z=h_actual)
    P = [40, 0, h_actual];
    W_punto = 0;
    
    % 4. Biot-Savart evaluado solo en el punto P
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
            %(B-S)
            V_ind = (G / (4 * pi)) * (c / mag_c2) * (dot(r1, r0)/norm(r1) - dot(r2, r0)/norm(r2));
            %Guardar downwash
            W_punto = W_punto + V_ind(3);
        end
    end
    
    %Teorema de Trefftz
    w_num = W_punto / 2;
    
    % Proporcionalidad marco T.
    L_num = L_inf * (w_inf / w_num);
end