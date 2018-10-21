# Leyes de escalamiento para calcular exponentes críticos

Siguiendo este [artículo](doi.org/10.1016/j.susc.2008.10.037), debemos construir muestras cuadradas de espín 1.0 tipo Ising. Con esto en mente, vamos a construir un script que nos construya las muestras para **vegas**. Después, debemos construir los respectivos archivos .json que contengan los parámetros. Luego, simulamos con **vegas** y por último, crearemos otro script para analizar los datos y obtener aquellos que se necesitan.

## Creación de script para construir muestras (build_sample.py)

Inicialmente, debemos definir la longitud (length) del sistema y la constante de intercambio (jex):

```python
length = 10
jex = 1.0
```

Luego creamos los sitios con *itertools*:

```python
sites = list(itertools.product(range(length), range(length)))

```

El siguiente paso consiste en general el diccionario de vecinos (nbhs):

```python
nbhs = {}
for i, site in enumerate(sites):
    x, y = site
    nbh1 = sites.index(((x + 1) % length, y))
    nbh2 = sites.index(((x - 1) % length, y))
    nbh3 = sites.index((x, (y + 1) % length))
    nbh4 = sites.index((x, (y - 1) % length))
    nbhs[i] = (nbh1, nbh2, nbh3, nbh4)

num_interactions = numpy.sum([len(nbhs[i]) for i in range(N)])
assert(num_interactions == 4*N)
```

Y, por último, exportamos el archivo de la muestra en el formato de **vegas**:

```python
with open("length_%i_.dat" % length, mode="w") as file:
    file.write("{} {} {}\n".format(N, num_interactions, 1))
    file.write("generic\n")

    for i in range(N):
        file.write("{} {} {} {} {} {} {} {} {} {}\n".format(
            i, *sites[i], 0.0, 1.0, 0, 0, 0, "generic", "flip"))

    for i in range(N):
        for j in nbhs[i]:
            file.write("{} {} {}\n".format(
                i, j, jex))


print("Sample with length = %i was created !!!" % length)
```


## Creación de script para construcción de archivos *json* (make_json.py)

Primero que todo, vamos a enlistar todos los archivos de las muestras para contruir los archivos *json* con base en ellos. Por otra parte, para el promedio estadístico, debemos hacer varias simulaciones, con diferentes números aleatorios, para cada muestra.

Para enlistar todos los archivos de muestra generados con el script anterior, vamos a hacer uso de la librería [glob](https://docs.python.org/3.5/library/glob.html):

```python
samples = sorted(glob("*.dat"))
```

Ahora, definimos el número de simulaciones independientes que se deben hacer por muestra:

```python
num_simulations = 5
```

Finalmente, hacemos un *for* anidado con el fin de construir diccionarios con los parámetros de las simulaciones y luego guardarlos en un archivo *json* por medio de la librería [json](https://docs.python.org/3.5/library/json.html). Es útil imprimir el comando que debe ser ingresado a la terminal para ejecutar **vegas**:

```python
for sample in samples:
    for i in range(num_simulations):
        parameters = dict(
            sample=sample,
            seed=numpy.random.randint(10000, 1000000),
            temperature=dict(start=5.0, final=0.01, points=200),
            mcs=10000,
            out=sample.replace(".dat", "sim_%i_.h5" % (i + 1)),
            kb=1.0
            )
        
        json_filename = sample.replace(".dat", "sim_%i_.json" % (i + 1))
        with open(json_filename, mode="w") as json_file:
            json.dump(parameters, json_file)

            print("time vegas %s > %s" % (json_filename, json_filename.replace(".json", ".log")))
```

## Creación de script para analizar los resultados (analyzer.py)

Vamos a crear un script para analizar los resultados y obtener las cantidad estadísticas de interés:
- <img src="https://latex.codecogs.com/png.latex?$\left<E\right>$" title="$\left<E\right>$" />: Energía media.
- <img src="https://latex.codecogs.com/png.latex?$C_{v}$" title="$C_{v}$" />: Calor específico.
- <img src="https://latex.codecogs.com/png.latex?$M$" title="$M$" />: Magnetización media.
- <img src="https://latex.codecogs.com/png.latex?$\chi$" title="$\chi$" />: Susceptibilidad magnética.
- <img src="https://latex.codecogs.com/png.latex?$V_{n}$" title="$V_{n}$" />: Cumulante de orden n del parámetro de orden.

Basándonos en este [artículo](doi.org/10.1016/j.susc.2008.10.037), estas cantidades son definidas como:

<img src="https://latex.codecogs.com/png.latex?$$\left<E\right>&space;=&space;\left<&space;\mathcal{H}&space;\right>$$" title="$$\left<E\right> = \left< \mathcal{H} \right>$$" />

<img src="https://latex.codecogs.com/png.latex?$$C_{v}&space;=&space;\frac{1}{k_B&space;T^2}&space;\left(&space;\left<E^2\right>&space;-&space;\left<E\right>^2&space;\right)$$" title="$$C_{v} = \frac{1}{k_B T^2} \left( \left<E^2\right> - \left<E\right>^2 \right)$$" />

<img src="https://latex.codecogs.com/png.latex?$$M&space;=&space;\left|\sum_i&space;\sigma_i&space;\right&space;|&space;$$" title="$$M = \left|\sum_i \sigma_i \right | $$" />

<img src="https://latex.codecogs.com/png.latex?$$\chi&space;=&space;\frac{1}{k_B&space;T}&space;\left(&space;\left<M^2\right>&space;-&space;\left<M\right>^2&space;\right)$$" title="$$\chi = \frac{1}{k_B T} \left( \left<M^2\right> - \left<M\right>^2 \right)$$" />

<img src="https://latex.codecogs.com/png.latex?$$V_{n}&space;=&space;\left<E\right>&space;-&space;\frac{\left<M^n&space;E\right>}{\left<M^n\right>}$$" title="$$V_{n} = \left<E\right> - \frac{\left<M^n E\right>}{\left<M^n\right>}$$" />