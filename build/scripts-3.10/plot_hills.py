#!python
import click
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergyLandscape


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot")
@click.option("--output", "-o", default=".", help="Output directory for figures")
@click.option("--time_resolution", "-tr", default=None, help="Number of decimal places for time values")
def main(file: str, output: str, time_resolution: int):

    hills = pl.read_as_pandas(file)
    landscape = FreeEnergyLandscape(hills)
    figures = landscape.get_hills_figures(time_resolution=time_resolution)

    for key, value in figures.items():
        save_dir = output + "/Walker_" + str(key) + ".pdf"
        value.write_image(save_dir)


if __name__ == "__main__":
    main()
