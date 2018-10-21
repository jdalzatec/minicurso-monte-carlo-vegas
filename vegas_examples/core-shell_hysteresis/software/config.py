import numpy

temps = []
fields = []
for T in numpy.linspace(1.0, 0.001, 100):
    temps.append(T)
    fields.append(1.0)

for H in numpy.linspace(1, -1, 100, endpoint=False):
    fields.append(H)
    temps.append(temps[-1])

for H in numpy.linspace(-1, 1, 101):
    fields.append(H)
    temps.append(temps[-1])


params = {
    "sample": "sample.dat",
    "anisotropy": "anisotropy.dat",
    "out": "results.h5",
    "kb": 1.0,
    "mcs": 1000,
    "seed": 696969,
    "field": fields,
    "temperature": temps
}

import json
json.dump(params, open("config.json", mode="w"))