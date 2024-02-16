import h5py

class IO:
    def __init__(self):
        pass

    def save_to_hdf5(self, grid, params, omega0s, omega1s, omega2s, output):
        filename = f'../data/spectral_funcs_dJU_{params.dJU}_UIB_{params.UIB:.2f}_Mu_{params.muU}_M_{grid.M}_N_{grid.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            # Save omega0s, omega1s, omega2s
            spectral_funcs_dset = f.create_dataset('spectral_funcs', data=self.spectral_funcs)
            f.create_dataset('omega0s', data=omega0s)
            f.create_dataset('omega1s', data=omega1s)
            f.create_dataset('omega2s', data=omega2s)
            f.create_dataset('output', data=output)

            # Create attributes for your parameters
            spectral_funcs_dset.attrs['UIB'] = self.UIB
            spectral_funcs_dset.attrs['Lx'] = self.grid.Lx
            spectral_funcs_dset.attrs['Ly'] = self.grid.Ly
            spectral_funcs_dset.attrs['muU'] = self.muU


    def save_to_hdf5_all(self, grid, params, dJU_values, en_values, omega0s, omega1s, omega2s, T11_values, T12_values, T21_values, T22_values, T22_SE_values, SE_SI):
        filename = f'../data/spectral_funcs_UIB_{params.UIB:.2f}_Mu_{params.muU}_M_{grid.M}_N_{grid.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            # Save omega0s, omega1s, omega2s
            for i, dJU in enumerate(dJU_values):
            # Create a group for this dJU value
                group = f.create_group(f'dJU_{dJU}')

                # Store the arrays in this group
                group.create_dataset('T11', data=T11_values[i, :])
                group.create_dataset('T12', data=T12_values[i, :])
                group.create_dataset('T21', data=T21_values[i, :])
                group.create_dataset('T22', data=T22_values[i, :])
                group.create_dataset('T22_SE', data=T22_SE_values[i, :])
                group.create_dataset('SE_SI', data=SE_SI)
                group.create_dataset('en', data=en_values)

            f.create_dataset('omega0s', data=omega0s)
            f.create_dataset('omega1s', data=omega1s)
            f.create_dataset('omega2s', data=omega2s)


            # Create attributes for your parameters
            f.attrs['UIB'] = params.UIB
            f.attrs['Lx'] = grid.Lx
            f.attrs['Ly'] = grid.Ly
            f.attrs['muU'] = params.muU



    def read_from_hdf5(self, UIB, muU, M, N):
        filename = f'../data/spectral_funcs_UIB_{UIB:.2f}_Mu_{muU}_M_{M}_N_{N}.hdf5'

        with h5py.File(filename, 'r') as f:
            # Read omega0s, omega1s, omega2s
            omega0s = f['omega0s'][:]
            omega1s = f['omega1s'][:]
            omega2s = f['omega2s'][:]

            # Read the data for each dJU value
            dJU_values = [key for key in f.keys() if key.startswith('dJU_')]
            data = {}
            for dJU in dJU_values:
                group = f[dJU]
                data[dJU] = {
                    'T11': group['T11'][:],
                    'T12': group['T12'][:],
                    'T21': group['T21'][:],
                    'T22': group['T22'][:],
                    'T22_SE': group['T22_SE'][:],
                    'SE_SI': group['SE_SI'][:],
                    'en': group['en'][:]
                }

            # Read the attributes
            UIB = f.attrs['UIB']
            Lx = f.attrs['Lx']
            Ly = f.attrs['Ly']
            muU = f.attrs['muU']

        return omega0s, omega1s, omega2s, data, UIB, Lx, Ly, muU
    
