clear; clc; close all;

% ============================================================
%  PANEL METHOD 2D: ESTUDIO DEL EFECTO SUELO
%  NACA 0012 | Método de Dirichlet con doble capa + wake panel
% ============================================================

%%
% =========================================================================
%  BLOQUE 1 — Datos
% =========================================================================

% -------------------------------------------------------------------------
%  Datos — Geometría NACA 0012
% -------------------------------------------------------------------------
N_cuerda = 60;
b        = 10.0;   % envergadura (solo para dimensionalizar)
c        = 1.0;   % cuerda
t_rel    = 0.12;

% -------------------------------------------------------------------------
%  Datos — Condiciones de vuelo
% -------------------------------------------------------------------------
Vmag   = 15;       % m/s
alpha  = 5;        % grados
rho    = 1.225;    % kg/m³
L_wake = 20*c;     % longitud del panel de estela (semi-infinita aproximada)


h_list = [0.10, 0.15, 0.20, 0.25, 0.30,0.35,0.40,0.45, 0.50,0.55, ... 
    0.60,0.65,0.70,0.75,0.80,0.85,0.90, 1.00,1.25, 1.50,1.75, 2.00];

%%
% =========================================================================
%  BLOQUE 2 — MAIN
% =========================================================================
[X,Z,N] = InicializarGeoNACA(N_cuerda,t_rel); % Geometría NACA

CL_vec = zeros(size(h_list));
L_vec  = zeros(size(h_list));
CM_vec = zeros(size(h_list));

for hh = 1:length(h_list)

    h = h_list(hh);

    % 1. Resolver potencial
    [mu, mu_w] = DirichletDoubletBEM(X,Z,N,h,Vmag,alpha,L_wake);

    % 2. Calcular fuerzas
    [CL,L,CM] = CL_Momentum_Fuerzas(X,Z,N,mu,Vmag,alpha,rho,h,b,c);

    % 3. Guardar resultados
    CL_vec(hh) = CL;
    L_vec(hh)  = L;
    CM_vec(hh) = CM;
end

%%
% =========================================================================
%  BLOQUE 3 — Gráficas
% =========================================================================

figure('Color','w','Position',[100 100 1300 420]);

subplot(1,3,1)
    plot(h_list, CL_vec, '-ob', 'LineWidth',2, 'MarkerFaceColor','b','Marker', 'none')
    hold on
    yline(CL_vec(end), '--k', sprintf('C_L^{free} = %.3f', CL_vec(end)), ...
          'LabelVerticalAlignment','bottom', 'LineWidth',1.2)
    grid on
    xlabel('h/c'); ylabel('C_L')
    title('Coeficiente de sustentación')
    xlim([0 h_list(end)+0.1])

subplot(1,3,2)
    plot(h_list, L_vec, '-sm', 'LineWidth',2, 'MarkerFaceColor','m','Marker', 'none')
    grid on
    xlabel('h/c'); ylabel('L [N]')
    title('Fuerza de sustentación')

subplot(1,3,3)
    plot(h_list, CM_vec, '-dg', 'LineWidth',2, 'MarkerFaceColor','g','Marker', 'none')
    grid on
    xlabel('h/c'); ylabel('C_M')
    title('Momento respecto a c/4')

sgtitle(sprintf('Panel Method 2D Corregido  —  \\alpha = %g°  (NACA 0012)', alpha), ...
        'FontSize',14)

%%
% =========================================================================
%  BLOQUE 4 — FUNCIONES LOCAES
% =========================================================================

function [X,Z,N] = InicializarGeoNACA(N_cuerda,t_rel)

    beta   = linspace(0, pi, floor(N_cuerda/2) + 1); %parmetrizac
    x_half = 0.5*(1 - cos(beta));
    y_half     = 5*t_rel*(0.2969*sqrt(x_half) - 0.1260*x_half - 0.3516*x_half.^2   + 0.2843*x_half.^3 - 0.1015*x_half.^4);

    X = [x_half(end:-1:1),  x_half(2:end)]; 
    Z = [y_half(end:-1:1),     -y_half(2:end)]; % Contorno cerrado

    % Orientación CCW %CCW = sentido contrahorario xd
    SArea = 0.5*sum(X(1:end-1).*Z(2:end) - X(2:end).*Z(1:end-1));
    if SArea < 0; X = fliplr(X); Z = fliplr(Z); end

    N = length(X) - 1;   % número de paneles
end

function [xc,tx,nx,zc,tz,nz,S,s_mid] = Prop_Paneles(X,Z,N)

    xc = zeros(N,1);  zc = zeros(N,1);
    tx = zeros(N,1);  tz = zeros(N,1);
    nx = zeros(N,1);  nz = zeros(N,1);
    S  = zeros(N,1);

    s_mid = zeros(N,1);
    arc   = 0;

    for i = 1:N
        dx    = X(i+1) - X(i);
        dz    = Z(i+1) - Z(i);

        S(i)  = sqrt(dx.^2 + dz.^2); %longitud del panel

        tx(i) = dx/S(i);
        tz(i) = dz/S(i); %vectores tangentes.

        % n = (tz, -tx)  →  vector normal.
        nx(i) =  tz(i);
        nz(i) = -tx(i);

        xc(i) = 0.5*(X(i) + X(i+1)); 
        zc(i) = 0.5*(Z(i) + Z(i+1)); %centro del panel

        % Longitud de arco acumulada hasta el punto de control de cada panel
        s_mid(i) = arc + 0.5*S(i);
        arc       = arc + S(i);
    end
end

function [mu, mu_w] = DirichletDoubletBEM(X,Z,N, h, Vmag,alpha,L_wake)

    [xc,~,~,zc,~,~,~,~] = Prop_Paneles(X,Z,N); % Propiedades geométricas de paneles
    
        % Posiciones de paneles a la altura h
        zc_h = zc + h;
        Zh   = Z(:) + h;   % nodos desplazados

        % Posición del borde de fuga (promedio entre nodo 1 y último nodo)
        xTE = 0.5*(X(1)  + X(end));
        zTE = 0.5*(Zh(1) + Zh(end)); %suavizar

        % Panel de estela: real e imagen
        % Normal deseada para estela real: (0,+1)  →  negamos dpPotential
        % Normal deseada para imagen estela: (0,-1) →  signo correcto por geometría
        wA  = [xTE,          zTE];    wB  = [xTE + L_wake,  zTE];
        wAi = [xTE,         -zTE];    wBi = [xTE + L_wake, -zTE];

    % ------------------------------------------------------------------
    %  Construcción de la matriz de influencia B  (N × N+1)
    %
    %  B(i,j) = potencial en pto. control i por panel j con μ=1
    %
    %  Reglas de signo (derivadas del método de imágenes):
    %    Panel real j    : +dpPotential(ri, Aj, Bj)
    %    Panel imagen j  : -dpPotential(ri, Aj_img, Bj_img)
    %
    %       [la imagen requiere negación: normal imagen = (nx,-nz),
    %        pero dpPotential usa right-of-tangent que da -(nx,-nz)]
    %
    %    Estela real     : -dpPotential(ri, wA, wB)
    %       [right-of-tangent da (0,-1); necesitamos (0,+1), negamos]
    %
    %    Imagen estela   : +dpPotential(ri, wAi, wBi)
    %       [right-of-tangent da (0,-1) = normal imagen correcta, OK]
    % ------------------------------------------------------------------

        B = zeros(N, N+1);

        for i = 1:N
            ri = [xc(i), zc_h(i)];

            for j = 1:N
                Aj  = [X(j),   Zh(j)  ];    Bj  = [X(j+1),   Zh(j+1)  ];
                Aji = [X(j),  -Zh(j)  ];    Bji = [X(j+1),  -Zh(j+1)  ];

                % Panel real (excluir auto-influencia; el término diagonal
                % 1/2 ya representa el salto en el límite interior)
                if i ~= j
                    B(i,j) = B(i,j) + dpPotential(ri, Aj, Bj);
                end

                % Panel imagen: siempre se incluye (está en z<0, lejos del cuerpo)
                B(i,j) = B(i,j) - dpPotential(ri, Aji, Bji);
            end

            % Columna de estela (real + imagen)
            B(i, N+1) = -dpPotential(ri, wA, wB) + dpPotential(ri, wAi, wBi);
        end

    % ------------------------------------------------------------------
    %  Sistema (N+1) × (N+1)
    %
    %  Filas 1..N  — BC Dirichlet: φ_int = 0
    %    (1/2)μ_i - Σ_j B(i,j)μ_j = φ∞(xi)
    %
    %  Fila N+1    — Condición de Kutta: μ_w = μ_1 - μ_N
    %    μ_1 - μ_N - μ_w = 0
    %
    %  El sistema resuelve simultáneamente los N coeficientes de doblete
    %  y la intensidad de la estela μ_w (circulación libre).
    % ------------------------------------------------------------------
    
        A_sys = zeros(N+1, N+1);
        for i = 1:N
            A_sys(i, :) = -B(i, :);   % -B fuera de diagonal
            A_sys(i, i) =  0.5;        % término de salto de la doble capa
        end

        % Kutta
        A_sys(N+1, 1)   =  1;    % μ_upper_TE  (panel 1)
        A_sys(N+1, N)   = -1;    % μ_lower_TE  (panel N)
        A_sys(N+1, N+1) = -1;    % μ_w

        % RHS: potencial de la corriente incidente φ∞ = V∞ · x
        phi_inf = Vmag*(xc*cosd(alpha) + zc_h*sind(alpha));
        RHS     = [phi_inf; 0];   % última fila: Kutta homogénea

        % Resolución
        mu_all = A_sys \ RHS;
        mu     = mu_all(1:N);
        mu_w = mu_all(N+1);   % disponible pero no se usa de momento
    
end

function phi = dpPotential(r_eval, r_A, r_B)

% =========================================================================
%    (ξ, η) son coordenadas locales del punto de campo:
%    ξ = proyección sobre tangente
%    η = proyección sobre normal (derecha del tangente)
%
%  Condición de salto verificada:
%    η → 0⁺ :  φ → +1/2   (exterior, lado de la normal)
%    η → 0⁻ :  φ → −1/2   (interior)
%    Salto  :  [φ] = 1 = μ  ✓
% =========================================================================

    dx = r_B(1) - r_A(1);
    dz = r_B(2) - r_A(2);
    L  = hypot(dx, dz);
    if L < 1e-14; phi = 0; return; end

    t_x = dx/L;   t_z = dz/L;

    % Coordenadas locales del punto de campo
    Rx = r_eval(1) - r_A(1);
    Rz = r_eval(2) - r_A(2);
    xi  =  Rx*t_x + Rz*t_z;        % coord. tangencial
    eta =  Rx*t_z - Rz*t_x;        % coord. normal  n=(t_z, -t_x)

    % Regularización: nunca se llama para i==j (auto-término), pero
    % paneles imagen cercanos pueden tener |η| pequeño con h/c pequeño.
    if abs(eta) < 1e-10
        eta = sign(eta + eps)*1e-10;
    end

    phi = (atan(xi/eta) - atan((xi - L)/eta)) / (2*pi);
end

function [Vt,Cp] = Vel_tangYCp(Vmag,mu,s_mid,N)
    % =========================================================================
    %  Debido a Dirichlet:
    %    φ_exterior_en_superficie = μ
    %
    %  Por lo tanto:
    %    V_t = ∂φ_ext/∂s = ∂μ/∂s
    %
    %  Se evalúa mediante diferencias finitas centrales sobre s_mid.
    % =========================================================================

    Vt = zeros(N,1);
    for i = 1:N
        if i == 1
            Vt(i) = (mu(2)  - mu(1)  ) / (s_mid(2)  - s_mid(1)  );
        elseif i == N
            Vt(i) = (mu(N)  - mu(N-1)) / (s_mid(N)  - s_mid(N-1));
        else
            Vt(i) = (mu(i+1) - mu(i-1)) / (s_mid(i+1) - s_mid(i-1));
        end
    end

    Cp = 1 - (Vt / Vmag).^2;
end

function [CL_h,L_h,CM_h]= CL_Momentum_Fuerzas(X,Z,N,mu,Vmag,alpha,rho,h,b,c)

    q_dyn  = 0.5*rho*Vmag^2;
    S_ref  = b*c;

    [xc,~,nx,zc,~,nz,S,s_mid] = Prop_Paneles(X,Z,N);
    [~ , Cp] = Vel_tangYCp(Vmag,mu,s_mid,N);
    zc_h = zc + h;

    Fp_prime = [0, 0];
    Mp_prime = 0;

    xref     = 0.25*c;
    zref     = h;

    for i = 1:N
        dF        = Cp(i)*q_dyn*S(i)*[nx(i), nz(i)];
        Fp_prime  = Fp_prime + dF;
        r_arm     = [xc(i)-xref, zc_h(i)-zref];
        Mp_prime  = Mp_prime + (r_arm(1)*dF(2) - r_arm(2)*dF(1));
    end

    F_total = Fp_prime * b;
    M_total = Mp_prime * b;

    eL           = [-sind(alpha), cosd(alpha)];
    L_val        = dot(F_total, eL);
    CL_h   = L_val / (q_dyn*S_ref);
    L_h    = L_val;
    CM_h   = M_total / (q_dyn*S_ref*c);

end

