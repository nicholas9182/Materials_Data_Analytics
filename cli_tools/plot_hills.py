#!/usr/bin/env python3
import click
import plumed as pl
from datetime import datetime
from analytics.metadynamics.free_energy import FreeEnergyLandscape


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot", type=str)
@click.option("--output", "-o", default=".", help="Output directory for figures", type=str)
@click.option("--time_resolution", "-tr", default=None, help="Number of decimal places for time values", type=int)
@click.option("--height_power", "-hp", default=1, help="Power to raise height of hills for easier visualisation", type=float)
def main(file: str, output: str, time_resolution: int, height_power: float):
    """
    cli tool to plot hill heights for all walkers, as well as the value of their CV
    :param file: the location of the HILLS file
    :param output: folder in which to put images
    :param time_resolution: how to bin the t axis for faster plotting
    :param height_power: power to raise hills too for easier visualisation
    :return: saved figures
    """
    hills = pl.read_as_pandas(file)
    landscape = FreeEnergyLandscape(hills)
    figures = landscape.get_hills_figures(time_resolution=time_resolution, height_power=height_power)

    for key, value in figures.items():
        key = str(key)
        save_dir = output + "/Walker_" + key + ".pdf"
        value.write_image(save_dir, scale=2)
        current_time = datetime.now().strftime("%H:%M:%S")
        click.echo(f"{current_time}: Made Walker_{key}.pdf in {output}", err=True)


if __name__ == "__main__":
    main()
