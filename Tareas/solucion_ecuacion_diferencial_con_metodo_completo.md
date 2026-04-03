# 📘 Resolución de Ecuación Diferencial

## 🔧 Problema
Resolver la ecuación diferencial:

\[
y'' + 3y' + \frac{25}{16}y = \cos x
\]

con condiciones iniciales:

\[
y(0) = 5, \quad y'(0) = 3
\]

---

# 🧠 1. Solución general

La solución general es:

\[
y(x) = y_h + y_p
\]

---

# ⚙️ 2. Solución homogénea

## Ecuación asociada

\[
y'' + 3y' + \frac{25}{16}y = 0
\]

## Ecuación característica

\[
\lambda^2 + 3\lambda + \frac{25}{16} = 0
\]

## Solución

\[
\lambda = -\frac{3}{2} \pm \frac{\sqrt{11}}{4}
\]

## Resultado

\[
y_h = C_1 e^{\left(-\frac{3}{2} + \frac{\sqrt{11}}{4}\right)x} + C_2 e^{\left(-\frac{3}{2} - \frac{\sqrt{11}}{4}\right)x}
\]

---

# 🔥 3. Solución particular (método complejo)

## Paso 1: reemplazo

\[
\cos x = \text{Re}(e^{ix})
\]

Resolvemos:

\[
y'' + 3y' + \frac{25}{16}y = e^{ix}
\]

## Paso 2: ansatz

\[
y_p = C e^{ix}
\]

## Paso 3: derivadas

\[
y_p' = iCe^{ix}, \quad y_p'' = -Ce^{ix}
\]

## Paso 4: sustitución

\[
C(-1 + 3i + \frac{25}{16})e^{ix} = e^{ix}
\]

\[
C\left(\frac{9}{16} + 3i\right) = 1
\]

## Paso 5: resolver C

\[
C = \frac{1}{\frac{9}{16} + 3i}
\]

Multiplicando por el conjugado:

\[
C = \frac{144}{2385} - \frac{768}{2385}i
\]

## Paso 6: parte real

\[
y_p = \frac{16}{265}\cos x + \frac{256}{795}\sin x
\]

---

# 📌 4. Solución general completa

\[
y(x) = C_1 e^{\lambda_1 x} + C_2 e^{\lambda_2 x} + \frac{16}{265}\cos x + \frac{256}{795}\sin x
\]

---

# 🧮 5. Aplicación de condiciones iniciales

## Condición 1: y(0)=5

\[
C_1 + C_2 + \frac{16}{265} = 5
\]

\[
C_1 + C_2 = \frac{1309}{265}
\]

---

## Condición 2: y'(0)=3

Derivando:

\[
y' = C_1\lambda_1 e^{\lambda_1 x} + C_2\lambda_2 e^{\lambda_2 x} - \frac{16}{265}\sin x + \frac{256}{795}\cos x
\]

Evaluando en 0:

\[
C_1\lambda_1 + C_2\lambda_2 + \frac{256}{795} = 3
\]

\[
C_1\lambda_1 + C_2\lambda_2 = \frac{2129}{795}
\]

---

# ⚡ 6. Cálculo de constantes

Sistema:

\[
C_1 + C_2 = S
\]

\[
C_1\lambda_1 + C_2\lambda_2 = T
\]

Método rápido:

\[
C_2 = S - C_1
\]

\[
C_1(\lambda_1 - \lambda_2) + S\lambda_2 = T
\]

\[
C_1 = \frac{T - S\lambda_2}{\lambda_1 - \lambda_2}
\]

\[
C_2 = S - C_1
\]

---

# ✅ 7. Resultado final

\[
y(x) = C_1 e^{\lambda_1 x} + C_2 e^{\lambda_2 x} + \frac{16}{265}\cos x + \frac{256}{795}\sin x
\]

con:

\[
\lambda_1 = -\frac{3}{2} + \frac{\sqrt{11}}{4}, \quad
\lambda_2 = -\frac{3}{2} - \frac{\sqrt{11}}{4}
\]

\[
C_1 = \frac{T - S\lambda_2}{\lambda_1 - \lambda_2}, \quad
C_2 = S - C_1
\]

---

# 🧠 Comentario final

La solución se compone de:

- Parte homogénea → comportamiento natural
- Parte particular → respuesta al forzamiento

El método complejo simplifica significativamente el cálculo de la solución particular.

