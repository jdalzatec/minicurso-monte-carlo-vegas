from glob import glob
import json
import numpy
import click

@click.command()
@click.option("--num_simulations", default=5, help="Amount of independent simulations for each sample.")
def main(num_simulations):

    samples = sorted(glob("*.dat"))

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
    

if __name__ == '__main__':
    main()