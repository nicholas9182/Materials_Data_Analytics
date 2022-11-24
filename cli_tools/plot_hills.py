import click
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergyLandscape


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot")
@click.option("--output", "-o", default=".", help="Output directory for figures")
def main(file, output):

    hills = pl.read_as_pandas(file)
    landscape = FreeEnergyLandscape(hills)
    figures = landscape.get_hills_figures()

    for key, value in figures.items():
        save_dir = output + "/Walker_" + str(key) + ".pdf"
        value.write_image(save_dir)


if __name__ == "__main__":
    main()
