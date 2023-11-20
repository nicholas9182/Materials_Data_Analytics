from analytics.metadynamics.free_energy import FreeEnergySpace

C12C8 = FreeEnergySpace.from_standard_directory("/Users/nsiemons/Drive/PostDoc/Projects/NDI/0005/005/")

data = (C12C8
        .get_reweighted_line_with_walker_error('DC1.min', bins=3, adaptive_bins=True)
        .get_data(with_metadata=True)
        )

data
