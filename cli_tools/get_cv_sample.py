#!/usr/bin/env python3
import click
import os
from Materials_Data_Analytics.metadynamics.free_energy import MetaTrajectory


@click.command()
@click.option("--colvar_file", "-f", type=str, help="COLVAR file to get sample from")
@click.option("--condition", "-c", type=str, multiple=True, help="condition, can take multiple conditions")
@click.option("--temperature", "-t", type=float, help="Temperature of the COLVAR", default=298)
@click.option("--sample_size", "-n", type=int, help="Number of samples to return", default=5)
@click.option("--traj_file", "-tr", type=str, help="Directory with the plumed.xtc", default=None)
@click.option("--out_group", "-g", type=str, help="group for which to output structures", default="non_Water")
@click.option("--output_structures", "-o", is_flag=True, default=False, help="output the structures?")
@click.option("--tpr_file", "-s", type=str, help="tpr file for the xtc trajectory")
@click.option("--ndx_file", "-nd", type=str, help="ndx file for the groups")
def main(colvar_file: str, condition: str, traj_file: str = None, sample_size: int = 5, temperature: float = 298,
         out_group: str = "non_Water", tpr_file: str = None, output_structures: bool = False, ndx_file: str = 'index.ndx'):
    """
    cli tool to get samples from trajectory under conditions
    :param colvar_file: the colvar file from which to get the sample
    :param condition: condition for which to get sample
    :param sample_size: how many samples to get
    :param temperature: temperature of the trajectory
    :param traj_file: directory with the plumed.xtc file
    :param out_group: group for which to output structures
    :param tpr_file: tpr file for outputting structures
    :param output_structures: output the structures?
    :param ndx_file: index file with groups
    :return:
    """
    data = MetaTrajectory(colvar_file=colvar_file, temperature=temperature).get_data()

    for c in condition:
        data = data.query(c)

    sample = data.sample(sample_size)

    if output_structures:
        counter = 1
        for i, r in sample.iterrows():
            time = r['time'] * 1000
            command = f"echo \"{out_group}\" | gmx trjconv -f {traj_file} -dump {time} -s {tpr_file} -pbc whole -o sample_{counter}.pdb -n {ndx_file}"
            os.system(command)
            counter += 1

    click.echo(sample)


if __name__ == "__main__":
    main()
