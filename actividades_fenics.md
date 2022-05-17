# Actividades con FEniCS. EDP y Métodos Numéricos

FEniCS es una biblioteca para la resolución de EDP mediante el método de los elementos finitos. Mantiene versiones C++ y Python. Su página web es <https://fenicsproject.org/>

FEniCS es muy ágil y está en continuo desarrollo. Aunque muy recientemente se ha publicado la versión FEniCSx, que será estándar a corto plazo, nos centraremos en las versiones anteriores (pues en nuestras pruebas se muestran más estables). Por tanto, estableceremos la siguiente documentación de referencia: <https://fenicsproject.org/olddocs/dolfin/latest/python/>

## Instalación

FEniCS está orientada a sistemas Unix. Su instalación en Windows requiere el uso de WSL (Windows subystem for Linux) o de contenedores como Docker. En sistemas Mac, puede utilizarse de forma estándar con Conda/Anaconda. En sistemas Linux/Ubuntu, simplemente hay que instalr el paquete de forma estándar (usando una *ppa*, un repositorio externo a Ubuntu).
- <https://fenicsproject.org/olddocs/dolfin/latest/python/installation.html>

- Como editor, se puede usar Atom o cualquier otro de los muchos editores adecuado para el lenguaje Python

## Primer ejemplo: Problema de Poisson con condiciones Dirichlet

- Como primer ejemplo, resolveremos la EDP de Poisson con condiciones de contorno Dirichlet:
$-\Delta u = f$  en  $\Omega$,
$u = 0$  sobre  $\Omega$.

- Podemos ver el código en el fichero <https://github.com/rrgalvan/edpymn/blob/main/src/fenics/0_poisson.py> para el caso en en el que $\Omega$ es el intervalo unidad, $(0,1)^2$ y $f=4$ (con lo que la solución exacta es $u(x,y)=x(x-1)y(y-1)$).
