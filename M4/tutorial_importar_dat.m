%% ========================================================================
%  TUTORIAL: IMPORTAR Y USAR ARCHIVOS .DAT DE PERFILES
%  Guía paso a paso con ejemplos prácticos
%% ========================================================================

clear; clc; close all;

fprintf('\n');
fprintf('================================================================\n');
fprintf('  TUTORIAL: IMPORTAR PERFILES DESDE ARCHIVOS .DAT\n');
fprintf('================================================================\n\n');

%% PASO 1: CREAR ARCHIVO DE EJEMPLO
fprintf('PASO 1: Creando archivo de ejemplo\n');
fprintf('-----------------------------------\n');

% Para este tutorial, vamos a crear archivos .dat de ejemplo
% En la práctica, descargarías estos de UIUC o Airfoil Tools

create_example_dat_files();

fprintf('\n✓ Archivos de ejemplo creados:\n');
fprintf('  - naca0012_example.dat\n');
fprintf('  - naca4412_example.dat\n');
fprintf('  - custom_airfoil.dat\n\n');

pause(1);

%% PASO 2: CARGAR UN ARCHIVO .DAT BÁSICO
fprintf('PASO 2: Cargando archivo .dat básico\n');
fprintf('-------------------------------------\n');

% Método más simple
[x, y, name] = AirfoilLoader.loadFromFile('naca0012_example.dat');

fprintf('✓ Archivo cargado exitosamente\n');
fprintf('  Nombre: %s\n', name);
fprintf('  Puntos: %d\n', length(x));
fprintf('  Rango X: [%.4f, %.4f]\n', min(x), max(x));
fprintf('  Rango Y: [%.4f, %.4f]\n\n', min(y), max(y));

pause(1);

%% PASO 3: VISUALIZAR LA GEOMETRÍA
fprintf('PASO 3: Visualizando la geometría\n');
fprintf('----------------------------------\n');

AirfoilLoader.visualize(x, y, name);

fprintf('✓ Gráfica generada\n');
fprintf('  - Borde de ataque (LE) marcado en verde\n');
fprintf('  - Borde de fuga (TE) marcado en rojo\n\n');

pause(2);

%% PASO 4: SEPARAR SUPERFICIES (Extradós/Intradós)
fprintf('PASO 4: Separando extradós e intradós\n');
fprintf('--------------------------------------\n');

[x_upper, y_upper, x_lower, y_lower] = AirfoilLoader.separateSurfaces(x, y);

% Visualizar separación
figure('Position', [100, 100, 800, 400]);
plot(x_upper, y_upper, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Extradós');
hold on;
plot(x_lower, y_lower, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Intradós');
plot(x_upper(1), y_upper(1), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
text(x_upper(1)+0.05, y_upper(1), 'TE', 'FontSize', 12);
xlabel('x/c', 'FontSize', 12);
ylabel('y/c', 'FontSize', 12);
title(sprintf('Superficies Separadas: %s', name), 'FontSize', 14);
legend('Location', 'best');
axis equal;
grid on;

fprintf('✓ Superficies separadas\n\n');

pause(2);

%% PASO 5: CALCULAR ESTADÍSTICAS DEL PERFIL
fprintf('PASO 5: Calculando estadísticas\n');
fprintf('--------------------------------\n');

stats = AirfoilLoader.computeStatistics(x, y);

fprintf('\n');

pause(1);

%% PASO 6: COMPARAR MÚLTIPLES PERFILES
fprintf('PASO 6: Comparando múltiples perfiles\n');
fprintf('--------------------------------------\n');

% Cargar varios perfiles
airfoil_files = {
    'naca0012_example.dat'
    'naca4412_example.dat'
    'custom_airfoil.dat'
};

figure('Position', [100, 100, 1200, 400]);

for i = 1:length(airfoil_files)
    [x_i, y_i, name_i] = AirfoilLoader.loadFromFile(airfoil_files{i});
    
    subplot(1, 3, i);
    plot(x_i, y_i, 'b-', 'LineWidth', 2);
    hold on;
    
    % Marcar puntos clave
    [~, idx_le] = min(x_i);
    plot(x_i(idx_le), y_i(idx_le), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot(x_i(1), y_i(1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    
    % Calcular espesor
    t_c = max(y_i) - min(y_i);
    
    xlabel('x/c');
    ylabel('y/c');
    title(sprintf('%s\nt/c ≈ %.1f%%', name_i, t_c*100), 'FontSize', 11);
    axis equal;
    grid on;
    xlim([-0.1, 1.1]);
end

sgtitle('Comparación de Perfiles Importados', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('✓ Comparación visual generada\n\n');

pause(2);

%% PASO 7: USAR CON EL SIMULADOR (Aproximación)
fprintf('PASO 7: Usando perfil importado con simulador\n');
fprintf('----------------------------------------------\n');

% 1. Cargar y Estadísticas
[x_naca, y_naca, name_naca] = AirfoilLoader.loadFromFile('naca4412_example.dat');
stats_naca = AirfoilLoader.computeStatistics(x_naca, y_naca);

t_c_approx = stats_naca.t_max_value;
f_c_approx = stats_naca.camber_max;

fprintf('\nParámetros extraídos de %s:\n', name_naca);
fprintf('  t/c estimado: %.3f\n', t_c_approx);
fprintf('  f/c estimado: %.3f\n', f_c_approx);

% 2. Crear Solver
h = 0.25;
alpha = 6;
V_inf = 50;

% USAMOS 'thickness' y 'camber' (no 'thickness_ratio')
solver = JoukowskiGroundEffect(h, alpha, V_inf, ...
    'thickness', t_c_approx, ...
    'camber', f_c_approx, ...
    'airfoil', name_naca); % Si este da error, borra esta línea

fprintf('\nResultados con efecto suelo:\n');
fprintf('  Altura (h/c): %.2f\n', h);
fprintf('  Ángulo (α):   %.1f°\n', alpha);
fprintf('  CL:           %.4f\n', solver.CL);

% 3. Visualizar (Ahora llamamos a la función de la clase que acabamos de arreglar)
solver.plotAirfoilWithGround();

pause(2);
%% PASO 8: EXPORTAR Y GUARDAR
fprintf('\nPASO 8: Exportar y guardar datos\n');
fprintf('--------------------------------\n');

% Guardar en formato MATLAB (.mat)
AirfoilLoader.saveAsMatlab(x_naca, y_naca, name_naca, 'my_airfoil.mat');

% También puedes exportar a CSV
T = table(x_naca, y_naca, 'VariableNames', {'x', 'y'});
writetable(T, 'airfoil_coords.csv');

fprintf('✓ Datos exportados a:\n');
fprintf('  - my_airfoil.mat\n');
fprintf('  - airfoil_coords.csv\n\n');

pause(1);

%% PASO 9: WORKFLOW COMPLETO DE PRINCIPIO A FIN
fprintf('PASO 9: Workflow completo ejemplo\n');
fprintf('----------------------------------\n');

fprintf('\nEJEMPLO: Analizar perfil descargado de UIUC\n\n');

% Simular descarga de archivo (en realidad usaremos el de ejemplo)
dat_file = 'naca4412_example.dat';

fprintf('1. Archivo descargado: %s\n', dat_file);

% Cargar
[x_wf, y_wf, name_wf] = AirfoilLoader.loadFromFile(dat_file);
fprintf('2. ✓ Geometría cargada (%d puntos)\n', length(x_wf));

% Analizar
stats_wf = AirfoilLoader.computeStatistics(x_wf, y_wf);
fprintf('3. ✓ Estadísticas calculadas\n');

% Visualizar
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(x_wf, y_wf, 'b-', 'LineWidth', 2);
xlabel('x/c'); ylabel('y/c');
title('Geometría Original');
axis equal; grid on;

% Crear solver aproximado
subplot(1,3,2);
solver_wf = JoukowskiGroundEffect(0.2, 5, 60, ...
    'thickness', stats_wf.t_max_value, ...
    'camber', stats_wf.camber_max);
% --- CÓDIGO CORREGIDO PARA EL SUBPLOT CENTRAL ---

% 1. Recoger las 4 variables separadas (Extradós e Intradós)
[xj_u, xj_l, yj_u, yj_l] = solver_wf.getAirfoilGeometry(100);

% 2. Plotear ambas partes
plot(xj_u, yj_u, 'r-', 'LineWidth', 2); 
hold on;
plot(xj_l, yj_l, 'r-', 'LineWidth', 2);

% 3. Ajustes estéticos (esto déjalo como estaba)
xlabel('x/c'); ylabel('y/c');
title('Aproximación Joukowski');
axis equal; grid on;

% Comparar CL a diferentes alturas
subplot(1,3,3);
h_range = linspace(0.1, 1.0, 20);
CL_range = zeros(size(h_range));

for i = 1:length(h_range)
    s_temp = JoukowskiGroundEffect(h_range(i), 5, 60, ...
        'thickness', stats_wf.t_max_value, ...
        'camber', stats_wf.camber_max);
    CL_range(i) = s_temp.CL;
end

plot(h_range, CL_range, 'g-', 'LineWidth', 2);
xlabel('h/c'); ylabel('C_L');
title('CL vs Altura');
grid on;

sgtitle(sprintf('Análisis Completo: %s', name_wf), 'FontSize', 14, 'FontWeight', 'bold');

fprintf('4. ✓ Análisis de efecto suelo completado\n');
fprintf('5. ✓ Gráficas generadas\n\n');

%% RESUMEN FINAL
fprintf('\n');
fprintf('================================================================\n');
fprintf('  RESUMEN: PASOS PARA USAR ARCHIVOS .DAT\n');
fprintf('================================================================\n\n');

fprintf('1. Descargar archivo .dat de UIUC o Airfoil Tools\n');
fprintf('2. Cargar: [x,y,name] = AirfoilLoader.loadFromFile(''archivo.dat'')\n');
fprintf('3. Visualizar: AirfoilLoader.visualize(x, y, name)\n');
fprintf('4. Extraer parámetros: stats = AirfoilLoader.computeStatistics(x, y)\n');
fprintf('5. Usar con simulador:\n');
fprintf('   solver = JoukowskiGroundEffect(h, alpha, V, ...\n');
fprintf('       ''thickness'', stats.t_max_value, ...\n');
fprintf('       ''camber'', stats.camber_max)\n');
fprintf('6. Analizar resultados: solver.CL, solver.plotAirfoilWithGround()\n\n');

fprintf('================================================================\n');
fprintf('  TUTORIAL COMPLETADO\n');
fprintf('================================================================\n');

%% ========================================================================
%  FUNCIONES AUXILIARES
%% ========================================================================

function create_example_dat_files()
    %CREATE_EXAMPLE_DAT_FILES Crea archivos .dat de ejemplo
    
    % NACA 0012
    create_naca_dat('naca0012_example.dat', 0, 0, 0.12);
    
    % NACA 4412
    create_naca_dat('naca4412_example.dat', 0.04, 0.4, 0.12);
    
    % Perfil personalizado
    create_custom_dat('custom_airfoil.dat');
end

function create_naca_dat(filename, m, p, t)
    %CREATE_NACA_DAT Crea archivo .dat de perfil NACA 4-dígitos
    
    % Generar coordenadas
    n_points = 50;
    x = cosspace(0, 1, n_points); % Distribución coseno
    
    % Espesor
    yt = 5*t * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + ...
                0.2843*x.^3 - 0.1015*x.^4);
    
    % Línea media (camber)
    if m > 0 && p > 0
        yc = zeros(size(x));
        for i = 1:length(x)
            if x(i) <= p
                yc(i) = m/p^2 * (2*p*x(i) - x(i)^2);
            else
                yc(i) = m/(1-p)^2 * ((1-2*p) + 2*p*x(i) - x(i)^2);
            end
        end
    else
        yc = zeros(size(x));
    end
    
    % Superficies
    x_upper = x;
    y_upper = yc + yt;
    
    x_lower = flipud(x(:));
    y_lower = flipud(yc(:)) - flipud(yt(:));
    
    % Escribir archivo
    fid = fopen(filename, 'w');
    
    % Encabezado
    m_int = round(m*100);
    p_int = round(p*10);
    t_int = round(t*100);
    fprintf(fid, 'NACA %d%d%02d (Example)\n', m_int, p_int, t_int);
    
    % Coordenadas (empezando por borde de fuga, extradós)
    for i = 1:length(x_upper)
        fprintf(fid, '%10.6f  %10.6f\n', x_upper(i), y_upper(i));
    end
    
    % Intradós
    for i = 1:length(x_lower)
        fprintf(fid, '%10.6f  %10.6f\n', x_lower(i), y_lower(i));
    end
    
    fclose(fid);
end

function create_custom_dat(filename)
    %CREATE_CUSTOM_DAT Crea perfil personalizado
    
    n = 40;
    theta = linspace(0, 2*pi, n);
    
    % Forma tipo Joukowski modificada
    r = 0.5 + 0.1*cos(theta);
    x = r .* cos(theta) + 0.5;
    y = r .* sin(theta) * 0.15;
    
    % Normalizar
    x = (x - min(x)) / (max(x) - min(x));
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'Custom Airfoil (Example)\n');
    for i = 1:length(x)
        fprintf(fid, '%10.6f  %10.6f\n', x(i), y(i));
    end
    fclose(fid);
end

function x = cosspace(x1, x2, n)
    %COSSPACE Distribución coseno para mejor resolución en LE
    beta = linspace(0, pi, n);
    x = x1 + (x2 - x1) * (1 - cos(beta)) / 2;
end
