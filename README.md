# Ground_Effect

# Análisis Aerodinámico Computacional: Estudio del Efecto Suelo

Este repositorio contiene dos implementaciones en MATLAB diseñadas para estudiar el **Efecto Suelo**. Se abordan dos enfoques complementarios: un análisis bidimensional (2D) mediante el Método de Paneles y un análisis tridimensional (3D) mediante el Método de Red de Vórtices (VLM - *Vortex Lattice Method*). 

Ambos códigos utilizan el **Método de las Imágenes** para simular la presencia del suelo (plano de simetría aerodinámico) a diferentes alturas relativas (`h/c`).

---

## M2. Método de Paneles 2D (`Panel_Method.m`)

Este script evalúa el comportamiento de un perfil alar **NACA 0012** operando cerca del suelo. 

**Características Destacadas:**
* **Formulación Matemática:** Emplea una distribución de singularidades tipo **doblete** de intensidad constante por panel y resuelve el problema aplicando la **Condición de Contorno de Dirichlet** (potencial interno perturbado igual a cero).
* **Estela (Wake):** Implementa la condición de Kutta de forma implícita mediante un panel de estela (vórtice) semi-infinito que fuerza a que las presiones en el borde de fuga se igualen.
* **Simulación del Suelo:** Refleja la geometría del perfil y su estela creando una "imagen" especular por debajo del suelo, invirtiendo la dirección de la normal para cumplir la condición de no penetración.

### 📋 Desglose de Funciones Locales

* `InicializarGeoNACA`: Genera las coordenadas (x, z) del perfil NACA 0012 utilizando una parametrización de semi-coseno para agrupar más paneles cerca del borde de ataque y de fuga.
* `Prop_Paneles`: Calcula las propiedades geométricas de la discretización: longitudes de panel, puntos de control (centros), vectores tangentes y vectores normales.
* `DirichletDoubletBEM`: **[Motor Principal]** Construye la matriz de coeficientes de influencia aerodinámica sumando las contribuciones de los paneles reales, los paneles imagen, la estela real y la estela imagen. Resuelve el sistema lineal para obtener las intensidades de los dobletes ($\mu$).
* `dpPotential`: Evalúa el potencial inducido por un panel de doblete unitario en cualquier punto del espacio basándose en coordenadas locales normales y tangenciales.
* `Vel_tangYCp`: Deriva numéricamente la intensidad de los dobletes respecto a la longitud de arco para encontrar la velocidad tangencial y, posteriormente, el coeficiente de presión ($C_p$).
* `CL_Momentum_Fuerzas`: Integra la distribución de presiones a lo largo de los paneles para calcular la sustentación total ($L$), el coeficiente de sustentación ($C_L$) y el momento de cabeceo respecto al cuarto de cuerda ($C_M$).

---

## 2. Vortex Lattice Method 3D (`VLM_ground_effect_1.m`)

Este script amplía el análisis a tres dimensiones, evaluando un ala de envergadura finita y un perfil alar **NACA 0012**. Permite visualizar cómo el efecto suelo no solo aumenta la sustentación, sino que altera drásticamente la resistencia inducida debido a la interacción tridimensional de los vórtices.

**Características Destacadas:**
* **Distribución en Malla:** Utiliza una malla con distribución de coseno tanto en la dirección de la cuerda como en la envergadura, lo que garantiza alta fidelidad geométrica en las puntas del ala (wingtips) y bordes.
* **Herraduras de Vórtices:** Modela el ala mediante vórtices de herradura (segmento ligado al cuarto de cuerda + dos filamentos libres semi-infinitos de estela).
* **Corrección de Oswald:** Incluye un factor de corrección empírico para alas de bajo alargamiento (Aspect Ratio), mejorando la precisión en el cálculo de la resistencia inducida ($C_{Di}$).

### 📋 Desglose de Funciones Locales

* `Oswald`: Calcula empíricamente el factor de eficiencia de Oswald basándose en la fórmula de Raymer, aplicando ajustes si el ala tiene un alargamiento bajo.
* `Malla_VLM`: Discretiza la geometría del ala, determinando las coordenadas de los vórtices ligados (a 1/4 del panel) y de los puntos de control (a 3/4 del panel) usando distribuciones de coseno.
* `Matriz_free`: Calcula la matriz de influencia aerodinámica para **vuelo libre** (sin suelo), evaluando el efecto que cada vórtice de herradura tiene sobre todos los puntos de control.
* `horseshoe_velocity`: Modela la velocidad inducida completa de una herradura sumando tres segmentos: el vórtice ligado, la pierna de estela izquierda y la pierna de estela derecha.
* `compute_circulation`: **[Núcleo Físico]** Aplica la Ley de Biot-Savart vectorialmente para calcular la velocidad inducida por un segmento rectilíneo de vórtice. Incluye un radio de núcleo `EPSILON` para evitar singularidades matemáticas.
* `compute_forces`: Aplica el Teorema de Kutta-Joukowski para integrar las fuerzas sustentadoras basándose en la circulación ($\Gamma$) de cada panel, obteniendo $C_L$, $C_{Di}$ y $C_M$.
* `Main_loop`: **[Motor Principal]** Itera a través de una lista de alturas relativas (`h/c`). En cada paso, calcula la influencia de los vórtices imagen invertidos, suma esta matriz a la de vuelo libre y resuelve el sistema lineal para obtener la nueva distribución de circulación.
* `compute_delta_cp`: Calcula el salto del coeficiente de presión diferencial ($\Delta C_p$) sobre la superficie del ala, útil para generar visualizaciones topográficas (mapas 3D de presiones).
