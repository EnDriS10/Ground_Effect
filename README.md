
# Ground_Effect

### ¿Cómo encontrar qué?
El repositorio está dividido en 3 carpetas (M1, M2, M3) cada una correspondiente a uno de los métodos. A su vez, cada carpeta se divide en dos subcarpetas, una llamada "Codigo base" dónde esta el motor de cada simulación y otra carpeta llamada "Graficas" dónde encontrar alteraciones del motor principal del método para construir las diferentes graficas a distintas alturas, envergaduras o sacar los errores. En la carpeta de gráficas de M1 también se incluyen las gráficas experimentales y de CFD.


####################################################################################################################################################################################################################################################################################################################

# Análisis Aerodinámico Computacional: Estudio del Efecto Suelo

Este repositorio contiene tres implementaciones en MATLAB diseñadas para estudiar el **Efecto Suelo**. Se abordan dos enfoques complementarios:un ala fina unidimensional en un espacio tridimensional (3D) con vórtices de Biot-Savart, un análisis bidimensional (2D) mediante el Método de Paneles y un análisis tridimensional (3D) mediante el Método de Red de Vórtices (VLM - *Vortex Lattice Method*). 

Todos los códigos utilizan el **Método de las Imágenes** para simular la presencia del suelo (plano de simetría aerodinámico) a diferentes alturas relativas (`h/c`).

---
## M1. Método de ala fina biot-savart + Extensión dinámica RK2 ('M1_RK2.m')

En este sript se comienza partiendo de la definición de una serie de parámetros que son: la velocidad incidente del viento, la densidad del aire, la longitud de la estela, la circulación, la altura y la envergadura. En primer lugar se obtiene la solución teórica para la sustentación y downwash a partir de las ecuaciones mostradas en el informe. Tras esto se construye una malla tridimensional con la función meshgrid para construir posteriormente un campo de velocidades del aire. En el siguiente bloque se establecen los filamentos tanto del ala real como de una simétrica respecto a Z=0 (correspondiente al método imagen que hará posible la simulación del efecto suelo, también se crea un vector gamma que se encargará de que la circulación en el ala real tenga el sentido contrario a la presente en el ala espejo. Para finalizar el bloque se aplica Biot-Savart en un bucle anidado para construir el campo de velocidades. La creación de este campo de velocidades pemitirá la visualización volumétrica del sistema en futuros bloques. Por su parte el Lift y el downwash numérico se calculan recurriendo a una función local que en lugar de construir el campo de velocidades completo genera solo el punto donde se toma la medida, esto perimite medir a la distancia deseada con un coste de computación menor, pero el funcionamiento es conceptualmente idéntico al de la representación 3D. Una vez se obtiene el Lift y downwash numerico, se compara con el analítico para determinar la discrepancia que sale mostrada en consola como porcentaje de error. Después de estop viene la sección encargada de la visualización, dónde a partir de la función stream3, se simula la trayectoria que seguiría una partícula en el campo de velocidades del aire, esto se representa en una figura y se crea el gráfico 3D.

En el bloque siguiente, se definen nuevas variables como la masa del ala, un coeficiente de rozamiento de orden 1 respecto a v, una altura inicial, la aceleración gravitacional, el tiempo de simulación y el paso de este. A partir de aquí se aplica el método de Runge-Kutta de segundo orden, calculando L(h) a cada paso y sumándole el resto de fuerzas, esto genera una trayectoria que se muestra por pantalla.

### 📋 Desglose de Funciones Locales
'calcular_L_numérico' : esta función aplica biot-savart a partir de unos parámetros como los iniciales, solo que tan solo en un punt, por lo que tan solo necesita un bucle simple. El retorno de la función es el únicamente el valor L_num.

Citas y referencia:

Este código es una entrega original y propia, sin embargo, como todo
trabajo académico tiene sus referencias. 

En primer lugar, el cálculo del lift viene inspirado en [17] (ver en pdf).
Además, para aprender el funcionamiento de funciones nativas de Matlab como
stream3, resolver bugs y optimizar el código (con sus tiempos de
compilación) en lugar del LLM nativo de Matlab (Copilot), se usó [Google]
Gemini. No se usó en ningun caso para sustituir el esfuerzo académico ni la originalidad,
sino como una herramienta para explorar funcionalidades y potenciar  el
trabajo.


## M2. Método de Paneles 2D ('Panel_Method.m')

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

## M3. Vortex Lattice Method 3D ('VLM_ground_effect_1.m')

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
