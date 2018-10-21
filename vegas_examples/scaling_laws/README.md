# Leyes de escalamiento para calcular exponentes críticos

Siguiendo este [artículo](doi.org/10.1016/j.susc.2008.10.037), debemos construir muestras cuadradas de espín 1.0 tipo Ising. Con esto en mente, vamos a construir un script que nos construya las muestras para **vegas**. Después, debemos construir los respectivos archivos .json que contengan los parámetros. Luego, simulamos con **vegas** y por último, crearemos otro script para analizar los datos y obtener aquellos que se necesitan.

## Creación de script para construir muestras (build_sample.py)

Inicialmente, debemos definir la longitud (length) del sistema y la constante de intercambio (Jex):

```python
length = 10
jex = 1.0
```