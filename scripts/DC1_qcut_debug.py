from analytics.metadynamics.free_energy import FreeEnergySpace

C20 = FreeEnergySpace.from_standard_directory("/Users/nsiemons/Drive/PostDoc/Projects/NDI/0005/008/")

data = (C20
        .get_reweighted_line_with_walker_error('CM3', bins=100, adaptive_bins=True)
        .get_data(with_metadata=True)
        )

print(data)
