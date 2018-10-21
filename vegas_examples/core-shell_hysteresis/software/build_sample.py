import numpy
from itertools import product
from matplotlib import pyplot
from collections import defaultdict

Rc = 7.5
R = 10.5

ss = sc = 1.0
jcc = 1.0
kc = 0.1 * jcc
ks = 0.5 * jcc
jint = -0.5 * jcc
jss = -0.5 * jcc
update_policy = "adaptive"


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

all_positions = numpy.array(sites)
core_positions = numpy.array(core_sites)
shell_positions = numpy.array(shell_sites)
core_positions = core_positions[core_positions[:, 1] == 0]
shell_positions = shell_positions[shell_positions[:, 1] == 0]


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

nhbs = defaultdict(list)
for site in sites:
    x, y, z = site
    for dx, dy, dz in [(1, 0, 0), (-1, 0, 0),
                       (0, 1, 0), (0, -1, 0),
                       (0, 0, 1), (0, 0, -1)]:
        nhb = ((x + dx), (y + dy), (z + dz))
        if nhb in sites:
            nhbs[site].append(nhb)

for site in sites:
    assert len(nhbs[site]) <= 6
    for nhb in nhbs[site]:
        assert numpy.linalg.norm(numpy.array(site) - numpy.array(nhb)) == 1.0
        assert site in nhbs[nhb]


types = dict()
for site in sites:
    prefix = "core" if site in core_sites else "shell"
    for nhb in nhbs[site]:
        nhb_prefix = "core" if nhb in core_sites else "shell"
        if prefix != nhb_prefix:
            prefix += "_interface"
            break
    types[site] = prefix


positions = defaultdict(list)
for site in sites:
    positions[types[site]].append(site)

print(positions.keys())

core = numpy.array(positions["core"])
shell = numpy.array(positions["shell"])
core_interface = numpy.array(positions["core_interface"])
shell_interface = numpy.array(positions["shell_interface"])
core = core[core[:, 1] == 0]
shell = shell[shell[:, 1] == 0]
core_interface = core_interface[core_interface[:, 1] == 0]
shell_interface = shell_interface[shell_interface[:, 1] == 0]

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


anisotropy_axis = dict()
field_axis = dict()
for site in sites:
    anisotropy_axis[site] = (0.0, 0.0, 1.0)
    field_axis[site] = (0.0, 0.0, 1.0)



num_interactions = 0
for site in sites:
    num_interactions += len(nhbs[site])
num_sites = len(sites)


sample_file = open("sample.dat", mode="w")
anisotropy_file = open("anisotropy.dat", mode="w")

sample_file.write("{} {} {}\n".format(num_sites, num_interactions, len(set(types.values()))))
print(num_sites, num_interactions, len(set(types.values())))


for t in sorted(set(types.values())):
    sample_file.write("{}\n".format(t))
    print(t)

for site in sites:
    i = sites.index(site)
    t = types[site]
    sample_file.write("{} {} {} {} {} {} {} {} {} {}\n".format(i, *site, spin[t], *field_axis[site], t, update_policy))
    anisotropy_file.write("{} {} {} {}\n".format(*anisotropy_axis[site], kan[t]))

for site in sites:
    t = types[site]
    for nhb in nhbs[site]:
        nhb_t = types[nhb]
        sample_file.write("{} {} {}\n".format(
            sites.index(site), sites.index(nhb),
            jex[(t.split("_")[0], nhb_t.split("_")[0])]))

sample_file.close()
anisotropy_file.close()