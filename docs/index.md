# Métodos Numéricos

## Método de Euler

El método de Euler se basa en la expansión de Taylor de la función $x(t)$. Tenemos
$$
\text{Expansión Taylor} \Rightarrow x(t+h) = x(t) + h\frac{dx}{dt} + \overbrace{ \frac{h^2}{2} \frac{d^2x}{dt^2} } ^{\epsilon} + O(h^3).
$$
Para avanzar en el tiempo la función por un paso $h$ lo suficientemente pequeño se utiliza la ecuación
$$
\boxed{x(t + h) = x(t) + hf(x,t).}
$$
El error asociado con la aproximación **está ligado a la cantidad de veces que hagamos la aproximación**. Lo podemos estimar con

$$
\sum\epsilon = \sum_{k=0}^{N-1}\frac{h^2}{2}\left. \frac{d^2x}{dt^2} \right|_{x_k, t_k} = \frac{h}{2}\sum_{k=0}^{N-1}h\left.\frac{df}{dt}\right|_{x_k, t_k}\\
\approx \frac{h}2\int_a^b\frac{df}{dt}d t = \frac{h}{2}\left[f_b - f_a\right].
$$

donde se asume que $N = (b-a)/h$.

* El algoritmo toma la siguiente forma:
  - Empezar con $t = t_0$, $x = x_0$
  - Discretizar el tiempo en pasos temporales de forma equidistante con espaciamiento $h$, donde cada punto en el tiempo está denotado con $t_i$
  - Para cada punto en el tiempo encontrar $x$ utilizando el resultado de la iteración previa: $x_i = x_{i-1} + hf(x_{i})$

## Método de Runge-Kutta 2$^{\rm do}$ Orden (RK2)

El método RK2 utiliza el punto medio para evaluar el método de Euler. Mientras que el método de Euler se aplica en el punto $t$ para evaluar la derivada para aproximar la función en el punto $x = t + h$, el método RK2 utiliza el punto medio $t + h/2$ alcanzando así una mejor aproximación para el mismo valor de $h$.

El método se deriva aplicando la serie de Taylor alrededor del punto medio $t + h/2$ para obtener el valor de la función en el punto $x(t + h)$. Se llega a:
$$
\boxed{x(t + h) = x(t) + hf[x(t + h/2), t + h/2] + O(h^3)}.
$$

El único problema es que requerimos conocer el valor de la función en el punto medio $x(t + h/2)$, el cual desconocemos.

Para aproximar este valor utilizamos el método de Euler con un paso $h/2$, $(x + h/2) = x(t) + \frac{h}{2}f(x,t)$. De esta manera, obtenemos las ecuaciones del método RK2:
* $k_{1} = hf(x,t),$
* $k_{2} = hf\left(x + \frac{k_{1}}{2},t + \frac{h}{2}\right)$
* $x(t + h) = x(t) + k_{2}$


## Método de Runge-Kutta de 4$^{\rm to}$ Orden (RK4)

La metodología anterior se puede aplicar aún a más puntos ubicados entre $x(t)$ y $x(t + h)$ realizando expansiones de Taylor. De esta forma se pueden agrupar términos de orden $h^3$, $h^4$, etc; para cancelar dichas expresiones. 

El problema de hacer esto es que las expresiones se vuelven más complicadas conforme incrementamos el orden de aproximación. En general, la regla de dedo es que el $4^{\rm to}$ orden corresponde al mejor compromiso entre complejidad y error de aproximación. Este método es el más utilizado comunmente para resolver ODEs. 

El álgebra para encontrar las ecuaciones de $4^{\rm to}$ orden es tediosa, pero el resultado final es
* $k_{1} = hf(x, t)$,
* $k_{2} = hf\left(x + \frac{k_{1}}{2}, t+\frac{h}2\right)$,
* $k_{3} = hf\left(x + \frac{k_{2}}{2}, t+\frac{h}2\right)$,
* $k_{4} = hf\left(x + k_{3}, t + h \right)$,
* $x(t+h) = x(t) + \frac{1}{6}(k_{1} + 2 k_{2} + 2k_{3} + k_{4})$.

Para la mayoría de aplicaciones, el método RK4 es el método de-facto para obtener soluciones. Es fácil de programar y devuelve resultados precisos. 

El error de aproximación es $O(h^5)$, mientras que el error global es aproximadamente del orden $O(h^4)$.

# Project Overview

::: ode

