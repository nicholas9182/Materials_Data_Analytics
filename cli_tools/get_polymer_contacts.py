#!/usr/bin/env python3
import click
import os
import pandas as pd
import MDAnalysis as mda
import networkx as nx
from MDAnalysis.analysis import contacts


@click.command()
@click.option("--reference_atoms", "-ref", type=str, help="Pymol atom selection language for contact group a")
@click.option("--selection_atoms", "-sel", type=str, help="Pymol atom selection language for contact group b")
@click.option("--tpr_file", "-s", type=str, default='Prod.tpr', help="tpr file for the MD universe")
@click.option("--xtc_file", "-f", type=str, default='Prod.xtc', help="trajectory file for the MD universe")
@click.option("--pol_num", "-n", type=int, default=100, help="Number of polymers in the universe")
@click.option("--output_tsv", "-o", default="contacts.tsv", help="Output file (should be a tsv)", type=str)
@click.option("--radius", "-r", default=4, help="radius in angstrom for the cut-off", type=float)
@click.option("--trajectory_slicer", "-ts", default=1, help="take every ts'th frame from the trajectory", type=int)
@click.option("--verbose", "-v", count=True)
@click.option("--sparse", "-sp", is_flag=True, default=False, help="Store as a sparse dataframe?")
def main(reference_atoms: str, selection_atoms: str, tpr_file: str, xtc_file: str, pol_num: int, output_tsv: str, radius: float,
         trajectory_slicer: int, verbose, sparse: bool = False):
    """
    Function to calculate edge data for a network using the polymers as a basis, and proximity as an edge
    :param reference_atoms: Pymol atom selection language for contact group a
    :param selection_atoms: Pymol atom selection language for contact group b
    :param tpr_file: tpr file for the MD universe
    :param xtc_file: trajectory file for the MD universe
    :param pol_num: Number of polymers in the universe
    :param output_tsv: Output file (should be a tsv)
    :param radius: radius for the cuttoff in angstrom
    :param trajectory_slicer: take every ts'th element from the time data in the trajectory
    :param verbose: how verbose to be
    :param sparse: output a sparse dataframe?
    :return: file with contact information
    """

    u = mda.Universe(tpr_file, xtc_file)
    connection_data = []
    pairs = []

    for r in range(0, pol_num):
        ref_group = u.atoms.fragments[r].select_atoms(reference_atoms)

        if verbose == 1:
            click.echo(f"Calculating contacts for polymer {r}", err=True)

        for s in range(0, pol_num):
            sel_group = u.atoms.fragments[s].select_atoms(selection_atoms)

            pairs.append((r, s))

            if (s, r) not in pairs and s != r:
                if verbose == 2:
                    click.echo(f"Calculating contacts between polymer {r} and polymer {s}", err=True)
                for ts in u.trajectory[::trajectory_slicer]:
                    if verbose == 3:
                        click.echo(f"Calculating contacts between polymer {r} and polymer {s} at time {u.trajectory.time}", err=True)
                    dist = contacts.distance_array(ref_group.positions, sel_group.positions)
                    n_contacts = contacts.contact_matrix(dist, radius).sum()

                    if sparse and n_contacts != 0:
                        connection_data.append(
                            pd.DataFrame({
                                'time': [u.trajectory.time],
                                'ref': [r],
                                'sel': [s],
                                'n': [n_contacts]
                            }))
                    elif not sparse:
                        connection_data.append(
                            pd.DataFrame({
                                'time': [u.trajectory.time],
                                'ref': [r],
                                'sel': [s],
                                'n': [n_contacts]
                            }))
                    else:
                        pass

    connection_data = pd.concat(connection_data)
    connection_data.to_csv(output_tsv, sep="\t")


if __name__ == "__main__":
    main()
