#!/usr/bin/env python3
import click
from datetime import datetime
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergySpace
from glob import glob


@click.command()
@click.option("--file", "-f", default="HILLS", help="Hills file to plot", type=str)
@click.option("--output", "-o", default="Figures/", help="Output directory for figures", type=str)
@click.option("--time_resolution", "-tr", default=6, help="Number of decimal places for time values", type=int)
@click.option("--height_power", "-hp", default=1, help="Power to raise height of _hills for easier visualisation", type=float)
@click.option("--bias_exchange", "-be", is_flag=True, default=False, help="Is this a bias-exchange simulation?")
def main(file: str, output: str, time_resolution: int, height_power: float, bias_exchange: bool = False):
    """
    cli tool to plot hill heights for all walkers, as well as the value of their CV. It also plots the average and max _hills deposited
    :param file: the location of the HILLS file
    :param output: folder in which to put images
    :param time_resolution: how to bin the t axis for faster plotting
    :param height_power: power to raise _hills too for easier visualisation
    :param bias_exchange: is this a bias-exchange simulation?
    :return: saved figures
    """
    if bias_exchange is True:
        file = [f for f in glob(file+"*")]
    landscape = FreeEnergySpace(file)
    figures = landscape.get_hills_figures(time_resolution=time_resolution, height_power=height_power)

    for key, value in figures.items():
        key = str(key)
        save_dir = output + "/Walker_" + key + ".pdf"
        value.update_traces(line_color='white')
        value.write_image(save_dir, scale=2)
        current_time = datetime.now().strftime("%H:%M:%S")
        click.echo(f"{current_time}: Made Walker_{key}.pdf in {output}", err=True)

    (landscape
     .get_average_hills_figure(time_resolution=time_resolution)
     .update_traces(line_color='white')
     .write_image(output + "/" + "hills_mean.pdf", scale=2)
     )
    current_time = datetime.now().strftime("%H:%M:%S")
    click.echo(f"{current_time}: Made hills_mean.pdf in {output}", err=True)

    (landscape
     .get_max_hills_figure(time_resolution=time_resolution)
     .update_traces(line_color='white')
     .write_image(output + "/" + "hills_max.pdf", scale=2)
     )
    current_time = datetime.now().strftime("%H:%M:%S")
    click.echo(f"{current_time}: Made hills_max.pdf in {output}", err=True)


if __name__ == "__main__":
    main()
