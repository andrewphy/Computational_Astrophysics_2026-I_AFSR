# Proyecto: Ecuación de Lane-Emden

Este proyecto resuelve numéricamente la ecuación de Lane-Emden para diferentes índices politrópicos usando el método de Runge-Kutta de cuarto orden (RK4).

## Descripción

La ecuación de Lane-Emden describe la estructura interna de una estrella bajo la aproximación politrópica:

$$
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n
$$

donde:
- $$\theta\$$ es la densidad adimensional
- $$\xi$$ es el radio adimensional
- $$n$$ es el índice politrópico
---
## Objetivo

El objetivo del proyecto es estudiar la estructura interna de estrellas bajo la aproximación politrópica, analizar cómo el índice politrópico \(n\) afecta las soluciones de la ecuación de Lane-Emden y comparar el caso \(n=3\) con el modelo solar estándar (Model S).

---
## Contexto físico de la ecuación de Lane-Emden

La ecuación de Lane-Emden surge del estudio del equilibrio hidrostático en una estrella, es decir, el balance entre dos efectos fundamentales:

La gravedad, que tiende a comprimir la materia hacia el centro.
El gradiente de presión, que actúa en sentido opuesto evitando el colapso.
### 1. Equilibrio hidrostático

La condición de equilibrio se expresa como:

$$
\frac{dP}{dr} = -\frac{G M(r)}{r^2} \rho
$$

Esta ecuación indica que la presión debe disminuir con el radio para contrarrestar la atracción gravitacional.

### 2. Conservación de masa

La masa contenida dentro de un radio (r) está dada por:

$$
\frac{dM}{dr} = 4\pi r^2 \rho
$$

Esta relación conecta la estructura de la estrella con su distribución de densidad.

### 3. Ecuación de estado politrópica

Para cerrar el sistema de ecuaciones, se asume una relación entre presión y densidad:

$$
P = K \rho^{1 + \frac{1}{n}}
$$

donde (n) es el índice politrópico. Esta aproximación simplifica la física del material estelar.

### 4. Cambio a variables adimensionales

Resolver directamente estas ecuaciones es complicado debido a la presencia de constantes físicas. Por ello, se introducen variables adimensionales:

$$
\rho = \rho_c \theta^n, \quad r = \alpha \xi
$$

Este cambio permite eliminar las constantes y obtener una ecuación universal.

### 5. Ecuación de Lane-Emden

Tras sustituir y simplificar, se obtiene:

$$
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n
$$

Esta ecuación describe la estructura interna de una estrella politrópica en términos de variables adimensionales.

### 6. Interpretación física

La solución $$(\theta(\xi))$$ representa cómo varía la densidad dentro de la estrella. El punto donde $$(\theta = 0)$$ define el radio de la estrella en unidades adimensionales.
## Ejemplo de salida

El siguiente gráfico muestra la comparación entre el modelo solar y el politropo \(n=3\):

![Comparación de densidad](https://github.com/andrewphy/Computational_Astrophysics_2026-I_AFSR/blob/main/Proyecto_Lane_Emden/Figuras/comparacion_densidad_n3_datosrealesvsmodelo.png)

## Características del código

- Método RK4 para integración numérica
- Manejo de la singularidad en \(\xi = 0\) mediante expansión de Taylor
- Determinación precisa de \(\xi_1\) usando interpolación cúbica (Cubic Spline)
- Verificación de derivadas mediante diferencias finitas


---

## Instalación

Clonar el repositorio
### Requerimientos

#### Manual
Usa "requirements.txt"

python3 -m venv entorno_ASTRO

source entorno_ASTRO/bin/activate

pip install -r requirements.txt

#### Automotico
Ejecutar el script de configuración:
chmod +x setup.sh
./setup.sh

## DISCLAIMER 
## Disclaimer

Los archivos PDF incluidos en este repositorio corresponden a material académico proporcionado por el docente del curso. La autoría y propiedad intelectual de dichos documentos pertenecen al profesor.

Estos materiales se incluyen únicamente con fines educativos y de referencia dentro del contexto de este proyecto. No se reclama autoría sobre ellos ni se pretende su redistribución con fines comerciales.

Los datos del archivo data_solar_real.csv fueron extraídos del modelo solar estándar Model S de Christensen-Dalsgaard et al. (1996), disponibles en: https://users-phys.au.dk/jcd/solar_models/
. Se utilizó el conjunto de variables que incluye perfiles de densidad y radio normalizado.
