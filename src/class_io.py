import h5py

class IO:
    def __init__(self, grid, params):
        self.spectral_funcs = params.spectral_funcs
        self.omegas = params.omegas
        self.inverse_kna_list = params.inverse_kna_list
        self.GS = params.GS
        self.Gamma = params.Gamma
        self.a_BB = params.a_BB
        self.NL = grid.NL
        self.NT = grid.NT
        self.Lambda = grid.Lambda
        self.T = grid.T


    def save_to_hdf5(self):
        if self.Gamma == []:
            filename = f'../data/spectral_funcs_a_BB_{self.a_BB:.2f}_NL_{self.NL}_NT_{self.NT}_Lambda_{self.Lambda}_T_{self.T}_No_Gamma.hdf5'
        else:
            filename = f'../data/spectral_funcs_a_BB_{self.a_BB:.2f}_NL_{self.NL}_NT_{self.NT}_Lambda_{self.Lambda}_T_{self.T}_With_Gamma.hdf5'
        with h5py.File(filename, 'w') as f:
            # Create datasets for your arrays
            spectral_funcs_dset = f.create_dataset('spectral_funcs', data=self.spectral_funcs)

            # Save the omegas
            omegas_dset = f.create_dataset('omegas', data=self.omegas)

            # Save the inverse_kna_list
            inverse_kna_list_dset = f.create_dataset('inverse_kna_list', data=self.inverse_kna_list)

            # save the GS
            GS_dset = f.create_dataset('GS', data=self.GS)

            # Create attributes for your parameters
            spectral_funcs_dset.attrs['Gamma'] = self.Gamma
            spectral_funcs_dset.attrs['a_BB'] = self.a_BB
            spectral_funcs_dset.attrs['NL'] = self.NL
            spectral_funcs_dset.attrs['NT'] = self.NT
            spectral_funcs_dset.attrs['Lambda'] = self.Lambda
            spectral_funcs_dset.attrs['T'] = self.T


    def read_from_hdf5(self):
        if self.Gamma == []:
            filename = f'../data/spectral_funcs_a_BB_{self.a_BB:.2f}_NL_{self.NL}_NT_{self.NT}_Lambda_{self.Lambda}_T_{self.T}_No_Gamma.hdf5'
        else:
            filename = f'../data/spectral_funcs_a_BB_{self.a_BB:.2f}_NL_{self.NL}_NT_{self.NT}_Lambda_{self.Lambda}_T_{self.T}_With_Gamma.hdf5'
        with h5py.File(filename, 'r') as f:
            # Read the spectral_funcs dataset
            self.spectral_funcs = f['spectral_funcs'][:]
            self.omegas = f['omegas'][:]
            self.inverse_kna_list = f['inverse_kna_list'][:]
            self.GS = f['GS'][:]

            # Read the attributes
            self.Gamma = f['spectral_funcs'].attrs['Gamma']
            self.a_BB = f['spectral_funcs'].attrs['a_BB']
            self.NL = f['spectral_funcs'].attrs['NL']
            self.NT = f['spectral_funcs'].attrs['NT']
            self.Lambda = f['spectral_funcs'].attrs['Lambda']
            self.T = f['spectral_funcs'].attrs['T']