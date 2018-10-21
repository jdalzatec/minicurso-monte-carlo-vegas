import numpy
import h5py
from matplotlib import pyplot
from collections import Counter

dataset = h5py.File("results.h5", mode="r")

print(list(dataset))

mcs = dataset.attrs.get("mcs")

types = [t.decode() for t in dataset.get("types")]
num_types = Counter(types)

print(num_types)

temperatures = dataset.get("temperature")[100:]
fields = dataset.get("field")[100:]

mag_z = dataset.get("magnetization_z")[100:, mcs//2:]
mag_z_by_type = {t: dataset.get("%s_z" % t)[100:, mcs//2:] for t in num_types}


mag_mean = numpy.mean(mag_z[:, ::10], axis=1) / numpy.sum(list(num_types.values()))
mag_mean_by_type = {t:  numpy.mean(
    mag_z_by_type[t][:, ::10], axis=1) / num_types[t] for t in num_types}


pyplot.figure(figsize=(8, 6))
pyplot.plot(fields, mag_mean, label=r"$M_{\rm total}$", lw=2)
for t, mag in mag_mean_by_type.items():
    pyplot.plot(fields, mag, label=r"$M_{\rm %s}$" % t.replace("_", "\ "), lw=2)

pyplot.xlabel(r"$H / J_{cc}$", fontsize=20)
pyplot.ylabel(r"$M$", fontsize=20)
pyplot.xlim(min(fields), max(fields))
pyplot.grid()
lgd1 = pyplot.legend(loc=9, fontsize=20,
              bbox_to_anchor=(0.0, 0.35, 1, 1), ncol=2, mode="expand")
pyplot.savefig("M_vs_H.png", bbox_extra_artists=(lgd1, ), bbox_inches='tight')
pyplot.close()


intercepts = list()
for i in range(len(fields) - 1):
    if numpy.sign(mag_mean[i]) != numpy.sign(mag_mean[i + 1]):
        intercepts.append((fields[i] + fields[i + 1]) * 0.5)
print(intercepts)

Hb = (intercepts[0] + intercepts[1]) * 0.5
Hc = numpy.abs(intercepts[0] - intercepts[1]) * 0.5
print(Hc, Hb)

pyplot.figure(figsize=(8, 6))
pyplot.plot(fields, mag_mean, label=r"$M_{\rm total}$", lw=2)
pyplot.axvline(Hb, color="crimson", lw=2, label=r"$H_{B}$")
pyplot.axvline(Hb - Hc, color="black", lw=2, label=r"$H_{B} - H_{C}$")
pyplot.axvline(Hb + Hc, color="orange", lw=2, label=r"$H_{B} + H_{C}$")
pyplot.xlabel(r"$H / J_{cc}$", fontsize=20)
pyplot.ylabel(r"$M$", fontsize=20)
pyplot.xlim(min(fields), max(fields))
pyplot.grid()
pyplot.legend(loc=4, fontsize=20)
pyplot.savefig("M_vs_H_lines.png")
pyplot.close()