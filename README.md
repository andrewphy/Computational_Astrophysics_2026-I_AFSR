# Proyecto: Ecuación de Lane-Emden

Este proyecto resuelve numéricamente la ecuación de Lane-Emden para diferentes índices politrópicos usando el método de Runge-Kutta de cuarto orden (RK4).

## Descripción

La ecuación de Lane-Emden describe la estructura interna de una estrella bajo la aproximación politrópica:

$$
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n
$$

donde:
- \(\theta\) es la densidad adimensional
- \(\xi\) es el radio adimensional
- \(n\) es el índice politrópico

El objetivo es encontrar la solución \(\theta(\xi)\) y el radio de la estrella \(\xi_1\) tal que:

$$
\theta(\xi_1) = 0
$$

---
## Objetivo

El objetivo del proyecto es estudiar la estructura interna de estrellas bajo la aproximación politrópica, analizar cómo el índice politrópico \(n\) afecta las soluciones de la ecuación de Lane-Emden y comparar el caso \(n=3\) con el modelo solar estándar (Model S).

---

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