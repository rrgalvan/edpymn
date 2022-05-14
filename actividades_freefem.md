# Actividades con FreeFEm++. EDP y Métodos Numéricos

FreeFEm++ un lenguaje y un motor para la resolución de EDP mediante el
MEF. Sitio web: \url{https://doc.freefem.org/}

## Instalación

- FreeFEm++ tiene licencia libre. Puede instalarse de la forma habitual en sistemas operativos Windows, Mac y GNU/Linux. Ver \url{https://doc.freefem.org/introduction/installation.html}. El lenguaje FreeFEm++ está inspirado en C/C++. Suele almacenarse en ficheros con la extensión `.edp` y es ejecutado por el intérprete FreeFEm++ (el programa FreeFEm++ instalado).

    - Además de FreeFEm++, se necesita instalar un editor que permita escribir programas en este lenguaje. Por ejemplo, se recomienda instalar \url{https://atom.io/}{Atom}, un editor de código abierto que está disponible para los sistemas operativos habituales. Vamos a ver cómo puede integrarse con FreeFEm++:

    - Una de las características de Atom es su capacidad para ser extendido mediante paquetes. Ver \href{https://flight-manual.atom.io/using-atom/sections/atom-packages/}{más información en su página web} (incluyendo el proceso para la Instalación de paquete). En particular, podemos instalar paquetes que le permitan (1) reconocer y colorear la sintaxis de FreeFEm++ y (2) utilizar FreeFEm++ para ejecutar el código que está siendo editado. En concreto:
      - `language-freefem-official`: paquete para la integración del
        lenguaje en Atom, en particular coloreado de la sintaxis
      - `freefem-runner`: para ejecutar ficheros (con la extensión `.epp`) desde el editor

    - Para ejecutar código FreeFEm++ usando `freefem-runner`. puede ser conveniente añadir la carpeta donde está situado el intérprete de FreeFEm++ a la variable `PATH` del sistema. Específicamente en sistemas Windows.
