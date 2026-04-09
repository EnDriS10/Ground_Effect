function error_numerico = M1(h, b)
   
    mostrar_graficos = true; 
    
    if mostrar_graficos
        figure('Color', 'k', 'Position', [100 100 1000 700]);
    end


    V_inf = 15;
    Gamma = 25;       
    L_estela = 80;   
    rho = 1.225;      
   -
    w_inf = Gamma / (-pi * b); 
  
    factor_h = (16 * (h/b)^2) / (1 + 16 * (h/b)^2);
    w_total_teorico = w_inf * factor_h;
    L_inf = rho * V_inf * Gamma * b;
    L_total_teorico = L_inf * (w_inf / w_total_teorico);
    
  
    
   
    
   
   
    L_numerico = calcular_L_numerico(h,b,Gamma,V_inf,rho,L_estela);
    error_numerico = 100 * abs((L_total_teorico - L_numerico) / L_total_teorico);
    
   



function L_num = calcular_L_numerico(h_actual, b, Gamma, V_inf, rho, L_estela)
   
    w_inf = Gamma / (-pi * b); 
    L_inf = rho * V_inf * Gamma * b;
    
   
    real_fil = {[L_estela, -b/2, h_actual; 0, -b/2, h_actual], ...
                [0, -b/2, h_actual; 0, b/2, h_actual], ...
                [0, b/2, h_actual; L_estela, b/2, h_actual]};
    espj_fil = {[L_estela, -b/2, -h_actual; 0, -b/2, -h_actual], ...
                [0, -b/2, -h_actual; 0, b/2, -h_actual], ...
                [0, b/2, -h_actual; L_estela, b/2, -h_actual]};
                
    filamentos = [real_fil, espj_fil];
    gammas = [Gamma, Gamma, Gamma, -Gamma, -Gamma, -Gamma];
    
   
    W_punto = 0; 
    
    
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
            V_ind = (G / (4 * pi)) * (c / mag_c2) * (dot(r1, r0)/norm(r1) - dot(r2, r0)/norm(r2));
            W_punto = W_punto + V_ind(3);
        end
    end
    
    
    w_num = W_punto / 2;
    
    
    L_num = L_inf * (w_inf / w_num);
end
end