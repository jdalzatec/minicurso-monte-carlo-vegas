from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import h5py
import click

@click.command()
# @click.argument("file")
@click.option("--file", default="length_10_sim_1_.h5")
def main(file):
    dataset = h5py.File(file, mode="r")
    mcs = dataset.attrs["mcs"]
    kb = dataset.attrs["mcs"]
    tau = mcs // 2
    mag = numpy.abs(dataset.get("magnetization_z")[:, tau:])

    energy = dataset.get("energy")[:, tau:]
    temperature = dataset.get("temperature")[:]
    N = len(dataset.get("types"))
    dataset.close()

    mean_energy = numpy.mean(energy, axis=1) / N
    specific_heat = numpy.std(energy, axis=1) ** 2 / (N * kb * temperature**2)
    mean_magnetization = numpy.mean(mag, axis=1) / N
    susceptibility = numpy.std(mag, axis=1) ** 2 / (N * kb * temperature)
    V1 = numpy.mean(energy, axis=1) - numpy.mean(mag**1 * energy, axis=1) / numpy.mean(mag**1, axis=1)
    V2 = numpy.mean(energy, axis=1) - numpy.mean(mag**2 * energy, axis=1) / numpy.mean(mag**2, axis=1)
    

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


    with open(file.replace(".h5", ".mean"), mode="w") as file_results:
        file_results.write("# T E Cv M Chi\n")
        for i, _ in enumerate(temperature):
            file_results.write("{} {} {} {} {} {} {}\n".format(
                temperature[i], mean_energy[i], specific_heat[i],
                mean_magnetization[i], susceptibility[i], V1[i], V2[i]
                ))

    print("File %s was analyzed !!!" % file)

if __name__ == '__main__':
    main()