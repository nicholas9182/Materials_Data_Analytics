#!/usr/bin/env python3
import click
from datetime import datetime
from analytics.metadynamics.free_energy import FreeEnergySpace


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot", type=str)
@click.option("--output", "-o", default="Figures/", help="Output directory for figures", type=str)
@click.option("--time_resolution", "-tr", default=6, help="Number of decimal places for time values", type=int)
@click.option("--height_power", "-hp", default=1, help="Power to raise height of hills for easier visualisation", type=float)
def main(file: str, output: str, time_resolution: int, height_power: float):
    """
    cli tool to plot hill heights for all walkers, as well as the value of their CV. It also plots the average and max hills deposited
    :param file: the location of the HILLS file
    :param output: folder in which to put images
    :param time_resolution: how to bin the t axis for faster plotting
    :param height_power: power to raise hills too for easier visualisation
    :return: saved figures
    """
    landscape = FreeEnergySpace(file)
    figures = landscape.get_hills_figures(time_resolution=time_resolution, height_power=height_power)

    for key, value in figures.items():
        key = str(key)
        save_dir = output + "/Walker_" + key + ".pdf"
        value.write_image(save_dir, scale=2)
        current_time = datetime.now().strftime("%H:%M:%S")
        click.echo(f"{current_time}: Made Walker_{key}.pdf in {output}", err=True)

    average_hills = landscape.get_average_hills_figure(time_resolution=time_resolution)
    average_hills.write_image(output + "/hills_mean.pdf", scale=2)
    current_time = datetime.now().strftime("%H:%M:%S")
    click.echo(f"{current_time}: Made hills_mean.pdf in {output}", err=True)
    max_hills = landscape.get_max_hills_figure(time_resolution=time_resolution)
    max_hills.write_image(output + "/hills_max.pdf", scale=2)
    current_time = datetime.now().strftime("%H:%M:%S")
    click.echo(f"{current_time}: Made hills_max.pdf in {output}", err=True)


if __name__ == "__main__":
    main()
