## larmojdsango

Este programa soluciona las ecuaciones de campo de Einstein acopladas mínimamente a 3 campos de materia: Klein-Gordon, Dirac y Proca. Este acople se hace a nivel de la acción de Einstein-Hilbert considerando los distintos lagrangianos.

El proceso de acople se puede ver en los siguientes artículos:
- Guzmán, F. S. (2009). The three dynamical fates of Boson Stars. Revista mexicana de física, 55(4), 321-326.
- Daka, E., Phan, N. N., & Kain, B. (2019). Perturbing the ground state of Dirac stars. Physical Review D, 100(8), 084042.

Para ocupar este progrma es necesario lo siguiente:
1. El archivo \'input.par\' guarda los parámetros de entrada, cada línea tiene especificada la variable. La variable \'campo0\' se refiere a uno de los parámetros libres del shootng, mientras que \'w\' se refiere al otro.
2. Una vez definido el parámetro de entrada se ejecuta la orden `make` sobre la carpeta raíz de este proyecto. Esto generará una carpeta llamada \'data\'.
3. En la carpeta \'data\' se encontrará una copia del archivo de parámetros \'input.par\' y un ejecutable llamado \'estrella\', se ejecuta este último con `./estrella` estando en la carpeta \'data\'.
4. Una vez terminado el proceso de ejecución se pueden graficar los distintos archivos generados.
