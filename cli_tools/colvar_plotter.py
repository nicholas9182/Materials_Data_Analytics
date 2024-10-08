#!/usr/bin/env python3
import click
import os
import pandas as pd
from visualisation.themes import custom_dark_template
import plotly.express as px
from datetime import datetime
from Materials_Data_Analytics.metadynamics.free_energy import MetaTrajectory, FreeEnergySpace


@click.command()
@click.option("--path", "-p", type=str, default='./', help="path with the files in it")
@click.option("--string_matcher", "-sm", type=str, default='COLVAR_REWEIGHT', help="string to match the files to be read in")
@click.option("--time_resolution", "-tr", default=3, help="Number of decimal places for time values", type=int)
@click.option("--output", "-o", default="Figures/", help="Output directory for figures", type=str)
def main(path: str, string_matcher: str, time_resolution: int, output: str):

    files = [path + f for f in os.listdir(path) if string_matcher in f and 'bck' not in f]
    space = FreeEnergySpace()

    click.echo("You are plotting the following files:")
    for f in files:
        click.echo(f)

    [space.add_metad_trajectory(MetaTrajectory(f)) for f in files]

    data = []
    for key, value in space.trajectories.items():
        new_data = value.get_data(with_metadata=True, time_resolution=time_resolution)
        data.append(new_data)
    data = pd.concat(data)

    if not space.opes:
        images = space.trajectories[0].cvs
    else:
        images = space.trajectories[0].cvs + ['zed', 'neff', 'nker']

    for image in images:
        figure = px.line(data,
                         x='time',
                         y=image,
                         template=custom_dark_template,
                         facet_row='walker',
                         labels={'time': 'Time [ns]', image: ''},
                         title=image,
                         width=1300, height=1000
                         )
        figure.update_traces(line={'width': 1}, line_color='white')
        figure.write_image(f"{output}/{image}.png")
        current_time = datetime.now().strftime("%H:%M:%S")
        click.echo(f"{current_time}: Made {image}.png in {output}", err=True)


if __name__ == "__main__":
    main()
