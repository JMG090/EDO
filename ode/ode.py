# ode/ode.py

"""Provee varios métodos numéricos.

Este módulo permite al usuario implementar métodos numéricos para resolver Ecuaciones Diferenciales Ordinarias

Examples: 
    >>> def fc(x,t): return -(x**3) + np.sin(t)
    >>> t1 = np.linspace(0.0, 5.0, 10)
    >>> euler(fc, t1, 0.0)
    array([ 0.        ,  0.        ,  0.29300855,  0.77691757,  1.06939524,
        0.83175805,  0.70976579,  0.40525197, -0.00931034, -0.54504162])
    >>> t2 = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    >>> rk2(fc, t2, 0.0)
    array([ 0.        ,  0.47942554,  0.87343374,  0.48707735,  0.01139971,
       -0.91669862])
    >>> rk4(fc, t1, 0.0)
    array([ 0.        ,  0.14998981,  0.52862388,  0.85362466,  0.92832171,
        0.84618259,  0.64771507,  0.32414365, -0.14797342, -0.64309438])


El módulo contiene las siguientes funciones:


- `euler(f, t, xi)` - Devuelve la solución de la EDO calculada con el método de Euler.
- `rk2(f, t, xi)` - Devuelve la solución de la EDO calculada con el método de Runge-Kutta de 2do orden.
- `rk4(f, t, xi)` - Devuelve la solución de la EDO calculada con el método de Runge-Kutta de 4to orden.

"""

def euler(f, t, xi):
    """Implementa el método de Euler.
    
    Examples:
        >>> t = np.linspace(0.0, 5.0, 10)
        >>> def fc(x,t): return -(x**3) + np.sin(t)
        >>> euler(fc, t, 0.0)
        array([ 0.        ,  0.        ,  0.29300855,  0.77691757,  1.06939524,
        0.83175805,  0.70976579,  0.40525197, -0.00931034, -0.54504162])


    Args:
       f (function): Ecuación diferencial ordinaria de la forma dx/dt = f(x,t)
       t (array): Grilla temporal de pasos equidistantes
       xi (int): Condición inicial para "t = 0"


    Returns:
       x (array): Devuelve un arreglo con las soluciones de la EDO en los tiempos t


    """

    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = xi

    for i in range(t.size-1):
        x[i+1] = x[i]+h*f(x[i],t[i])
    
    return x


def rk2(f, t, xi):
    """Implementa el método de Runge-Kutta de 2do orden.
    
    Examples:
        >>> t = np.linspace(0.0, 5.0, 10)
        >>> def fc(x,t): return -(x**3) + np.sin(t)
        >>> rk2(fc, t, 0.0)
        array([ 0.        ,  0.15234405,  0.54887121,  0.85911636,  0.88561732,
        0.79451062,  0.60406615,  0.28737967, -0.1879154 , -0.69147835])                    

    Args:
       f (function): Ecuación diferencial ordinaria de la forma dx/dt = f(x,t)
       t (array): Grilla temporal de pasos equidistantes
       xi (int): Condición inicial para t=0


    Returns:
       x (array): Devuelve un arreglo con las soluciones de la EDO en los tiempos t


    """


    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = xi

    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i]+(k1/2),t[i]+(h/2))
        x[i+1] = x[i] + k2
    
    return x


def rk4(f, t, xi):
    """Implementa el método de Runge-Kutta de 4to orden.

    Examples:
        >>> t = np.linspace(0.0, 5.0, 10)
        >>> def fc(x,t): return -(x**3) + np.sin(t)
        >>> rk4(fc, t, 0.0)
        array([ 0.        ,  0.14998981,  0.52862388,  0.85362466,  0.92832171,
        0.84618259,  0.64771507,  0.32414365, -0.14797342, -0.64309438])

    Args:
       f (function): Ecuación diferencial ordinaria de la forma dx/dt = f(x,t)
       t (array): Grilla temporal de pasos equidistantes
       xi (int): Condición inicial para t=0


    Returns:
       x (array): Devuelve un arreglo con las soluciones de la EDO en los tiempos t


    """


    h = t[1]-t[0]
    x = np.zeros(t.size)
    x[0] = xi

    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i]+(k1/2),t[i]+(h/2))
        k3 = h * f(x[i]+(k2/2),t[i]+(h/2))
        k4 = h * f(x[i]+k3,t[i]+h)
        x[i+1] = x[i]+(1/6)*(k1+2*k2+2*k3+k4)
    
    return x
