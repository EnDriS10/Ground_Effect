%% ========================================================================
%  EXTENSIÓN: PERFILES AERODINÁMICOS EN EFECTO SUELO
%  Transformación de Joukowski + Schottky-Klein
%% ========================================================================

classdef JoukowskiGroundEffect < GroundEffectAnalytical
    properties
        thickness_ratio  
        camber_ratio     
        airfoil_type     
        
        a_joukowski      % Ahora será el RADIO REAL del círculo (R)
        b_joukowski      % Centro del círculo
        
        thickness_correction  
        camber_correction     
    end
    
    methods
        function obj = JoukowskiGroundEffect(h, alpha, V_inf, varargin)
            obj = obj@GroundEffectAnalytical(h, alpha, V_inf);
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'thickness', 0.12);
            addParameter(p, 'camber', 0.0);
            addParameter(p, 'airfoil', 'Custom'); 
            addParameter(p, 'chord', 1.0);
            parse(p, varargin{:});
            
            % --- SALVAVIDAS: Proteger contra datos NaN/vacíos del AirfoilLoader ---
            t_val = p.Results.thickness;
            if isempty(t_val) || isnan(t_val), t_val = 0.12; end
            obj.thickness_ratio = t_val;
            
            c_val = p.Results.camber;
            if isempty(c_val) || isnan(c_val), c_val = 0.0; end
            obj.camber_ratio = c_val;
            
            ch_val = p.Results.chord;
            if isempty(ch_val) || isnan(ch_val), ch_val = 1.0; end
            if isprop(obj, 'chord')
                obj.chord = ch_val;
            end
            
            % --- MATEMÁTICAS CORREGIDAS ---
            % La constante de mapeo (c) debe ser un cuarto de la cuerda
            c_mapping = ch_val / 4; 
            
            real_shift = -0.8 * obj.thickness_ratio * c_mapping;
            imag_shift =  2.0 * obj.camber_ratio * c_mapping;
            
            obj.b_joukowski = complex(real_shift, imag_shift);
            
            % El radio real (R) debe calcularse usando Pitágoras para cerrar el perfil
            % R = sqrt((c - real(b))^2 + imag(b)^2)
            obj.a_joukowski = sqrt((c_mapping - real_shift)^2 + imag_shift^2);
            
            obj.Gamma = pi * obj.V_inf * ch_val * sin(obj.alpha);
        end
        
        function obj = computeJoukowskiParameters(obj)
            obj.a_joukowski = 1 + obj.thickness_ratio/4;
            obj.b_joukowski = obj.camber_ratio * obj.chord;
        end
        
        function obj = solveWithAirfoilCorrections(obj)
            obj = solve(obj);
            K_thickness = obj.computeThicknessCorrection();
            K_camber = obj.computeCamberCorrection();
            obj.thickness_correction = K_thickness;
            obj.camber_correction = K_camber;
            CL_flat_plate = obj.CL;
            obj.CL = CL_flat_plate * K_thickness * K_camber;
            obj.Gamma = obj.CL * 0.5 * obj.rho * obj.V_inf^2 * obj.chord / (obj.rho * obj.V_inf);
        end
        
        function K_t = computeThicknessCorrection(obj)
            t_c = obj.thickness_ratio;
            beta = 0.5; 
            K_t = 1 - beta * (t_c^2) / (1 + obj.h/obj.chord);
            ground_factor = exp(-obj.h/obj.chord);
            K_t = 1 - (1 - K_t) * (1 + ground_factor);
        end
        
        function K_c = computeCamberCorrection(obj)
            f_c = obj.camber_ratio;
            if f_c == 0, K_c = 1.0; return; end
            alpha_eff_camber = 2 * f_c; 
            ground_amplification = 1 + 1/(4*obj.h/obj.chord + 1);
            Delta_CL_camber = 2*pi*alpha_eff_camber * ground_amplification;
            K_c = 1 + Delta_CL_camber / obj.CL;
        end
        
        function [x_upper, x_lower, y_upper, y_lower] = getAirfoilGeometry(obj, n_points)
            c_mapping = obj.chord / 4;
            R = obj.a_joukowski;
            b = obj.b_joukowski;
            
            theta = linspace(0, 2*pi, n_points);
            zeta = R * exp(1i*theta) + b;
            z = zeta + c_mapping^2 ./ zeta; 
            
            x = real(z);
            y = imag(z);
            
            [~, idx_split] = min(x); 
            
            x_upper_raw = x(1:idx_split);
            y_upper_raw = y(1:idx_split);
            
            x_lower_raw = x(idx_split:end);
            y_lower_raw = y(idx_split:end);
            
            x_min = min(x);
            x_range = max(x) - x_min;
            if x_range == 0, x_range = 1; end 
            
            c = obj.chord;
            x_upper = (x_upper_raw - x_min) / x_range * c;
            y_upper = y_upper_raw / x_range * c;
            
            x_lower = (x_lower_raw - x_min) / x_range * c;
            y_lower = y_lower_raw / x_range * c;
        end
        
        function plotAirfoilWithGround(obj)
            figure('Name', 'Resultados Internos de la Clase (Joukowski)', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 500]);
            
            % --- 1. GEOMETRÍA ---
            subplot(1,3,1);
            [x_u, x_l, y_u, y_l] = obj.getAirfoilGeometry(200);
            
            % Extradós en Cian, Intradós en Amarillo, Suelo en Blanco Punteado
            plot(x_u, y_u + obj.h, 'c-', 'LineWidth', 2.5); hold on;
            plot(x_l, y_l + obj.h, 'y-', 'LineWidth', 2.5);
            line([-0.5, 1.5], [0, 0], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
            
            title('Geometría en Efecto Suelo', 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('x/c', 'FontSize', 12, 'FontWeight', 'bold'); 
            ylabel('y/c', 'FontSize', 12, 'FontWeight', 'bold');
            axis equal; grid on;
            
            max_y = max(y_u + obj.h);
            ylim([-0.2, max(max_y + 0.5, 0.0)]); 
            
            % --- 2. DISTRIBUCIÓN DE PRESIÓN (Cp) ---
            subplot(1,3,2);
            
            c_mapping = obj.chord / 4;
            R = obj.a_joukowski;
            b = obj.b_joukowski;
            
            theta = linspace(0, 2*pi, 200);
            z_circle = R * exp(1i*theta); 
            zeta = z_circle + b; 
            z_airfoil = zeta + c_mapping^2 ./ zeta; 
            
            term1 = obj.V_inf * exp(-1i*obj.alpha);
            term2 = -obj.V_inf * (R^2 ./ z_circle.^2) * exp(1i*obj.alpha);
            term3 = -1i * obj.Gamma / (2*pi) ./ z_circle;
            dW_dzeta = term1 + term2 + term3;
            
            dz_dzeta = 1 - (c_mapping ./ zeta).^2;
            
            V_local = zeros(size(theta));
            for i = 1:length(theta)
                if abs(dz_dzeta(i)) < 1e-4
                    V_local(i) = 0; 
                else
                    V_local(i) = abs(dW_dzeta(i) / dz_dzeta(i));
                end
            end
            
            Cp = 1 - (V_local / obj.V_inf).^2;
            
            % Línea principal del Cp en Cian brillante
            plot(real(z_airfoil)/obj.chord, Cp, 'c-', 'LineWidth', 2.5);
            
            title('Distribución de Presión (Cp)', 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('x/c', 'FontSize', 12, 'FontWeight', 'bold'); 
            ylabel('C_p', 'FontSize', 12, 'FontWeight', 'bold');
            set(gca, 'YDir', 'reverse'); 
            ylim([-5, 1.5]); 
            grid on;
            
            % --- 3. CL vs ALTURA ---
            subplot(1,3,3);
            h_vals = linspace(0.1, 1.0, 50); 
            cl_vals = zeros(size(h_vals));
            
            Gamma_orig = obj.Gamma;
            for i = 1:length(h_vals)
                factor_suelo = (obj.chord / (4 * h_vals(i)))^2;
                Gamma_temp = Gamma_orig * (1 + factor_suelo);
                cl_vals(i) = 2 * Gamma_temp / (obj.V_inf * obj.chord);
            end
            
            % Línea en Verde brillante, punto rojo más grande
            plot(h_vals, cl_vals, 'g-', 'LineWidth', 2.5); hold on;
            plot(obj.h, obj.CL, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
            
            title('Curva de Efecto Suelo', 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('Altura (h/c)', 'FontSize', 12, 'FontWeight', 'bold'); 
            ylabel('Coeficiente C_L', 'FontSize', 12, 'FontWeight', 'bold');
            grid on;
        end 
    end
end