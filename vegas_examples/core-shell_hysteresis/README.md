# Ciclos de histéresis para nanopartículas core/shell

Vamos a construir una nanopartícula core/shell con estructura cúbica simple usando un script de python.

Primero que todo, debemos incluir las librerías necesarias:

```python
import numpy
from itertools import product
from matplotlib import pyplot
from collections import defaultdict
```

Ahora, definimos los valores del radio del core <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;R_c\right)" title="\left ( R_c\right)" /> y de la nanopartícula <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;R\right)" title="\left ( R\right)" /> medidos en celdas unitarias.

```
Rc = 7.5
R = 10.5
```

Definimos el valor del espín para el core <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;S_c\right)" title="\left ( S_c\right)" /> y para el shell <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;S_s\right)" title="\left ( S_s\right)" />. Así mismo, debemos definir el valor de las constantes de anisotropía para el core <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;k_c\right)" title="\left ( k_c\right)" /> y para el shell <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;k_s\right)" title="\left ( k_s\right)" /> y los valores de las constantes de intercambio entre core-core <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;j_{cc}\right)" title="\left ( j_{cc}\right)" />, core-shell <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;j_{int}\right)" title="\left ( j_{int}\right)" /> y shell-shell <img src="https://latex.codecogs.com/gif.latex?\left&space;(&space;j_{ss}\right)" title="\left ( j_{ss}\right)" />. Por otra parte, debemos definir la [política de actualización](https://pcm-ca.github.io/vegas/spin-update-policies/) para la dirección de los momentos magnéticos.

```python
ss = sc = 1.0
jcc = 1.0
kc = 0.1 * jcc
ks = 0.5 * jcc
jint = -0.5 * jcc
jss = -0.5 * jcc
update_policy = "adaptive"
```

Creamos diccionarios para fácilmente devolver los parámetros anteriormente definidos:

```python
spin = {
    "core": sc,
    "core_interface": sc,
    "shell": ss,
    "shell_interface": ss,
       }

kan = {
    "core": kc,
    "core_interface": kc,
    "shell": ks,
    "shell_interface": ks,
      }

jex = {
    ("core", "core"): jcc,
    ("core", "shell"): jint,
    ("shell", "core"): jint,
    ("shell", "shell"): jss
      }

```

Ahora, creamos la lista de sitios y los identificamos como iones del core o del shell de acuerdo a <img src="https://latex.codecogs.com/gif.latex?R_c" title="R_c" />:

```python
sites = list()
core_sites = list()
shell_sites = list()
for site in product(range(-int(numpy.ceil(R+1)), int(numpy.ceil(R+1))),
                    range(-int(numpy.ceil(R+1)), int(numpy.ceil(R+1))),
                    range(-int(numpy.ceil(R+1)), int(numpy.ceil(R+1)))):
    dist = numpy.linalg.norm(site)
    if dist <= R:
        sites.append(site)
        if dist <= Rc:
            core_sites.append(site)
        else:
            shell_sites.append(site)

```

Convertimos las listas anteriores a arreglos de ```numpy``` para producir un corte transversal para <img src="https://latex.codecogs.com/gif.latex?y&space;=&space;0" title="y = 0" />:

```python
all_positions = numpy.array(sites)
core_positions = numpy.array(core_sites)
shell_positions = numpy.array(shell_sites)
core_positions = core_positions[core_positions[:, 1] == 0]
shell_positions = shell_positions[shell_positions[:, 1] == 0]
```

Generamos una gráfica del sección transversal usando diferente colores para los iones del core y del shell:

```python
pyplot.figure(figsize=(8, 8))
pyplot.scatter(core_positions[:, 0], core_positions[:, 2], s=100, color="silver", edgecolor="black")
pyplot.scatter(shell_positions[:, 0], shell_positions[:, 2], s=150, color="gray", edgecolor="black")
pyplot.grid()
pyplot.xlabel(r"$x$", fontsize=20)
pyplot.ylabel(r"$z$", fontsize=20)
pyplot.xlim(-R-1, R+1)
pyplot.ylim(-R-1, R+1)
pyplot.gca().set_aspect("equal")
pyplot.show()

```
![cross-section](https://pcm-ca.github.io/vegas/tutorials/system-building/building-a-core-shell-nanoparticle/output_16_0.png)

Identificamos los vecinos de cada sitio y los almacenamos en un diccionario:

```python
nhbs = defaultdict(list)
for site in sites:
    x, y, z = site
    for dx, dy, dz in [(1, 0, 0), (-1, 0, 0),
                       (0, 1, 0), (0, -1, 0),
                       (0, 0, 1), (0, 0, -1)]:
        nhb = ((x + dx), (y + dy), (z + dz))
        if nhb in sites:
            nhbs[site].append(nhb)
```

Podemos realizar algunas verificaciones. Cada sitio debe tener como máximo 6 vecinos y debe estar a una distancia de 1.0 de cada uno de ellos:

```python
for site in sites:
    assert len(nhbs[site]) <= 6
    for nhb in nhbs[site]:
        assert numpy.linalg.norm(numpy.array(site) - numpy.array(nhb)) == 1.0
        assert site in nhbs[nhb]

```

Ahora, creamos un diccionario para identifica el tipo de cada sitio, el cual puede ser ```core```, ```shell```, ```core_interface``` o ```shell_interface```, los cuales corresponde a los sitios localizados en el core, el shell, la interfaz del core y la interfaz del shell, respectivamente:

```python
types = dict()
for site in sites:
    prefix = "core" if site in core_sites else "shell"
    for nhb in nhbs[site]:
        nhb_prefix = "core" if nhb in core_sites else "shell"
        if prefix != nhb_prefix:
            prefix += "_interface"
            break
    types[site] = prefix

```

Creamos un diccionario para almacenas todos los sitios de cada tipo:

```python
positions = defaultdict(list)
for site in sites:
    positions[types[site]].append(site)
```

Verificamos que sólo hay cuatro tipo de iones:

```python
positions.keys()
```

```python
dict_keys(['core_interface', 'core', 'shell', 'shell_interface'])
```

De nuevo, convertimos la lista de sitios a arreglos de nu py para poder generar un gráfico de un corte transversal para <img src="https://latex.codecogs.com/gif.latex?y&space;=&space;0" title="y = 0" />:

```python
core = numpy.array(positions["core"])
shell = numpy.array(positions["shell"])
core_interface = numpy.array(positions["core_interface"])
shell_interface = numpy.array(positions["shell_interface"])
core = core[core[:, 1] == 0]
shell = shell[shell[:, 1] == 0]
core_interface = core_interface[core_interface[:, 1] == 0]
shell_interface = shell_interface[shell_interface[:, 1] == 0]
```

Generamos un gráfico de la sección transversal usando diferentes colores para las diferentes regiones:

```python
pyplot.figure(figsize=(8, 8))
pyplot.scatter(core[:, 0], core[:, 2], s=100, color="silver", edgecolor="black")
pyplot.scatter(shell[:, 0], shell[:, 2], s=150, color="gray", edgecolor="black")
pyplot.scatter(core_interface[:, 0], core_interface[:, 2], s=100, color="gold", edgecolor="black")
pyplot.scatter(shell_interface[:, 0], shell_interface[:, 2], s=150, color="red", edgecolor="black")
pyplot.grid()
pyplot.xlabel(r"$x$", fontsize=20)
pyplot.ylabel(r"$z$", fontsize=20)
pyplot.xlim(-R-1, R+1)
pyplot.ylim(-R-1, R+1)
pyplot.gca().set_aspect("equal")
pyplot.show()

```

![cross-section](https://pcm-ca.github.io/vegas/tutorials/system-building/building-a-core-shell-nanoparticle/output_30_0.png)

Definimos la anisotropía y el campo magnético externo para cada sitio, los cuales van a ser iguales a la dirección z <img src="https://latex.codecogs.com/gif.latex?\left(0,&space;0,&space;1&space;\right&space;)" title="\left(0, 0, 1 \right )" /> para todos los iones:

```python
anisotropy_axis = dict()
field_axis = dict()
for site in sites:
    anisotropy_axis[site] = (0.0, 0.0, 1.0)
    field_axis[site] = (0.0, 0.0, 1.0)
```

Contamos el número de interacciones, la cual es igual a la suma de la cantidad de vecinos de cada sitio, y el número de iones, el cual es la longitud de la lista de sitios:

```python
num_interactions = 0
for site in sites:
    num_interactions += len(nhbs[site])
num_sites = len(sites)
```

Creamos los archivos para almacenar las propiedades estructurales (samples.dat) y la anisotropía (anisotropy.dat):

```python
sample_file = open("sample.dat", mode="w")
anisotropy_file = open("anisotropy.dat", mode="w")
```

Escribimos en la primera línea de **sample_file** el número de sitios, interacciones y tipos:

```python
sample_file.write("{} {} {}\n".format(num_sites, num_interactions, len(set(types.values()))))
print(num_sites, num_interactions, len(set(types.values())))
```

```
4945 27576 4
```

Escribimos el tipo de los iones, uno en una línea diferente:

```python
for t in sorted(set(types.values())):
    sample_file.write("{}\n".format(t))
    print(t)
```

```
core
core_interface
shell
shell_interface
```

Escribimos los parámetros de cada sitio de acuerdo al [formato](https://pcm-ca.github.io/vegas/system-building/) establecido por **Vegas**:

```python
for site in sites:
    i = sites.index(site)
    t = types[site]
    sample_file.write("{} {} {} {} {} {} {} {} {} {}\n".format(i, *site, spin[t], *field_axis[site], t, update_policy))
    anisotropy_file.write("{} {} {} {}\n".format(*anisotropy_axis[site], kan[t]))
```

Escribimos las interacciones de intercambio entre cada par de vecinos:

```python
for site in sites:
    t = types[site]
    for nhb in nhbs[site]:
        nhb_t = types[nhb]
        sample_file.write("{} {} {}\n".format(
            sites.index(site), sites.index(nhb),
            jex[(t.split("_")[0], nhb_t.split("_")[0])]))
```

Cerramos los archivos:

```python
sample_file.close()
anisotropy_file.close()
```