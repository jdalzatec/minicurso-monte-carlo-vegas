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

<img src="https://latex.codecogs.com/png.latex?$$M&space;=&space;\frac{1}{N}\left|\sum_i&space;\sigma_i&space;\right&space;|&space;$$" title="$$M = \frac{1}{N}\left|\sum_i \sigma_i \right | $$" />

<img src="https://latex.codecogs.com/png.latex?$$\chi&space;=&space;\frac{1}{k_B&space;T}&space;\left(&space;\left<M^2\right>&space;-&space;\left<M\right>^2&space;\right)$$" title="$$\chi = \frac{1}{k_B T} \left( \left<M^2\right> - \left<M\right>^2 \right)$$" />

<img src="https://latex.codecogs.com/png.latex?$$V_{n}&space;=&space;\left<E\right>&space;-&space;\frac{\left<M^n&space;E\right>}{\left<M^n\right>}$$" title="$$V_{n} = \left<E\right> - \frac{\left<M^n E\right>}{\left<M^n\right>}$$" />

Primero que todo, debemo importar algunas librerías para graficación, manejo de arreglos numéricos y manejo de archivos HDF5

```python
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import h5py
```

Ahora, cargamos el archivo HDF5 donde están los resultados obnetidos por **vegas**. Supongamos que el archivo se llama *length_10_sim_1_.h5*, entonces:

```python
dataset = h5py.File(file, mode="r")
```

La variable *dataset* ahora contiene toda la información. Podemos obtener algunos parámetros desde el diccionario de atributos:

```python
mcs = dataset.attrs["mcs"]
kb = dataset.attrs["mcs"]
```

y definimos la variable **tau** como la mitad de los pasos Monte Carlo con el fin de despreciar los primeros **tau** pasos Monte Carlo para relajación:

```python
tau = mcs // 2
```

Además, podemos obtener la evolución de la magnetización y la energía. En este caso, vamos a cargar la magentización en z puesto que es un modelo de Ising con una política de actualización tipo *flip*:

```python
mag = numpy.abs(dataset.get("magnetization_z")[:, tau:])
energy = dataset.get("energy")[:, tau:]
```

Así mismo, cargamos el arreglo de temperaturas:

```python
temperature = dataset.get("temperature")[:]
```
También podemos calcular al número de iones como la longitud del arreglo de posiciones:

```python
N = len(dataset.get("positions"))
```

Con todo esto cargado, podemos cerrar el **dataset** puesto que no necesitamos cargar más datos:

```python
dataset.close()
```

Empleando la definición de las variables arriba mencionadas, calculamos los observables:

```python
mean_energy = numpy.mean(energy, axis=1)
specific_heat = numpy.std(energy, axis=1) ** 2 / (kb * temperature**2)
mean_magnetization = numpy.mean(mag, axis=1) / N
susceptibility = numpy.std(mag, axis=1) ** 2 / (kb * temperature)
V1 = numpy.mean(energy, axis=1) - numpy.mean(mag**1 * energy, axis=1) / numpy.mean(mag**1, axis=1)
V2 = numpy.mean(energy, axis=1) - numpy.mean(mag**2 * energy, axis=1) / numpy.mean(mag**2, axis=1)

```

Podemos realizar gráficas de los observables calculados:

```python
pdf = PdfPages(file.replace(".h5", ".pdf"))

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.plot(temperature, mean_energy, "-s", color="black")
ax2.plot(temperature, specific_heat, "-o", color="crimson")
ax.grid()
ax.set_xlabel(r"$T$", fontsize=20)
ax.set_ylabel(r"$E$", fontsize=20, color="black")
ax2.set_ylabel(r"$C_{v}$", fontsize=20, color="crimson")
pyplot.tight_layout()
pyplot.savefig(pdf, format="pdf")
pyplot.close()

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.plot(temperature, mean_magnetization, "-s", color="black")
ax2.plot(temperature, susceptibility, "-o", color="crimson")
ax.grid()
ax.set_xlabel(r"$T$", fontsize=20)
ax.set_ylabel(r"$M$", fontsize=20, color="black")
ax2.set_ylabel(r"$\chi$", fontsize=20, color="crimson")
pyplot.tight_layout()
pyplot.savefig(pdf, format="pdf")
pyplot.close()

fig = pyplot.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.plot(temperature, V1, "-s", color="black")
ax2.plot(temperature, V2, "-o", color="crimson")
ax.grid()
ax.set_xlabel(r"$T$", fontsize=20)
ax.set_ylabel(r"$V_1$", fontsize=20, color="black")
ax2.set_ylabel(r"$V_2$", fontsize=20, color="crimson")
pyplot.tight_layout()
pyplot.savefig(pdf, format="pdf")
pyplot.close()

pdf.close()
```

Y finalmente, guardamos los observables en un archivo, tal que puedan ser analizados posteriormente y aplicarles las leyes de escalamiento:

```python
with open(file.replace(".h5", ".mean"), mode="w") as file_results:
    file_results.write("# T E Cv M Chi\n")
    for i, _ in enumerate(temperature):
        file_results.write("{} {} {} {} {} {} {}\n".format(
            temperature[i], mean_energy[i], specific_heat[i],
            mean_magnetization[i], susceptibility[i], V1[i], V2[i]
            ))

print("File %s was analyzed !!!" % file)
```

## Aplicando leyes de escalamiento (scaling.py)