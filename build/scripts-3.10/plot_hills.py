#!python
import click
import plumed as pl
from datetime import datetime
from analytics.metadynamics.free_energy import FreeEnergyLandscape


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot", type=str)
@click.option("--output", "-o", default=".", help="Output directory for figures", type=str)
@click.option("--time_resolution", "-tr", default=None, help="Number of decimal places for time values", type=int)
def main(file: str, output: str, time_resolution: int):

    hills = pl.read_as_pandas(file)
    landscape = FreeEnergyLandscape(hills)
    figures = landscape.get_hills_figures(time_resolution=time_resolution)

    for key, value in figures.items():
        key = str(key)
        save_dir = output + "/Walker_" + key + ".pdf"
        value.write_image(save_dir)
        current_time = datetime.now().strftime("%H:%M:%S")
        click.echo(f"{current_time}: Made Walker_{key}.pdf in {output}", err=True)


if __name__ == "__main__":
    main()
