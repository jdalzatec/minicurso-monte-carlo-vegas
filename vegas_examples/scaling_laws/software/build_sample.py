import numpy
from matplotlib import pyplot
import itertools
import click

@click.command()
@click.option("--length", "-l", default=50, help="Length of the system.")
@click.option("--jex", default=1.0, help="Exchange constant.")
def main(length, jex):
    N = length * length

    sites = []
    spins = []

    sites = list(itertools.product(range(length), range(length)))

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

if __name__ == '__main__':
    main()