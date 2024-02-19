import h5py
import numpy as np

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


    def save_to_hdf5_all_fixed_chemical_potentials(self, grid, params, dJU_values, en_values, omega0s, omega1s, omega2s, T11_values, T12_values, T21_values, T22_values, T22_SE_values, SE_SI):
        filename = f'../data/spectral_funcs_UIB_{params.UIB:.2f}_Mu_{params.muU}_M_{grid.M}_N_{params.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            f.create_dataset('en_values', data=en_values)
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
                group.create_dataset('SE_SI', data=SE_SI[i, :])

            f.create_dataset('omega0s', data=omega0s)
            f.create_dataset('omega1s', data=omega1s)
            f.create_dataset('omega2s', data=omega2s)
            f.create_dataset('dJU_values', data=dJU_values)

            # Create attributes for your parameters
            f.attrs['UIB'] = params.UIB
            f.attrs['Lx'] = grid.Lx
            f.attrs['Ly'] = grid.Ly
            f.attrs['muU'] = params.muU
            f.attrs['N'] = params.N

    def save_to_hdf5_all_fixed_density(self, grid, params, n0, dJU_values, muU_values, en_values, omega0s, omega1s, omega2s, T11_values, T12_values, T21_values, T22_values, T22_SE_values, SE_SI):
        filename = f'../data/spectral_funcs_UIB_{params.UIB:.2f}_n0_{n0}_M_{grid.M}_N_{params.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            f.create_dataset('en_values', data=en_values)
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
                group.create_dataset('SE_SI', data=SE_SI[i, :])

            f.create_dataset('omega0s', data=omega0s)
            f.create_dataset('omega1s', data=omega1s)
            f.create_dataset('omega2s', data=omega2s)
            f.create_dataset('dJU_values', data=dJU_values)
            f.create_dataset('muU_values', data=muU_values)

            # Create attributes for your parameters
            f.attrs['UIB'] = params.UIB
            f.attrs['Lx'] = grid.Lx
            f.attrs['Ly'] = grid.Ly
            f.attrs['n0'] = n0
            f.attrs['N'] = params.N

    def save_to_hdf5_perturbative(self, grid, params, dJU_values, omega0s, omega1s, omega2s):
        filename = f'../data/Perturbative_UIB_{params.UIB:.2f}_Mu_{params.muU}_M_{grid.M}_N_{params.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            # Save omega0s, omega1s, omega2s
            f.create_dataset('omega0s', data=omega0s)
            f.create_dataset('omega1s', data=omega1s)
            f.create_dataset('omega2s', data=omega2s)
            f.create_dataset('dJU_values', data=dJU_values)


            # Create attributes for your parameters
            f.attrs['UIB'] = params.UIB
            f.attrs['Lx'] = grid.Lx
            f.attrs['Ly'] = grid.Ly
            f.attrs['N'] = params.N

    def save_to_hdf5_fixed_density_qcorr(self, grid, params, desired_n0, dJU_values, mu_qcorr):
        filename = f'../data/fixed_density_line_chemical_potentials_n0_{desired_n0}_UIB_{params.UIB:.2f}_M_{grid.M}_N_{params.N}.hdf5'

        with h5py.File(filename, 'w') as f:
            f.create_dataset('mu_qcorr', data=mu_qcorr)
            f.create_dataset('dJU_values', data=dJU_values)

            # Create attributes for your parameters
            f.attrs['UIB'] = params.UIB
            f.attrs['Lx'] = grid.Lx
            f.attrs['Ly'] = grid.Ly
            f.attrs['muU'] = params.muU
            f.attrs['N'] = params.N
            f.attrs['n0'] = desired_n0
    
    def read_from_hdf5_fixed_density_qcorr(self, grid, params, desired_n0):
        filename = f'../data/fixed_density_line_chemical_potentials_n0_{desired_n0}_UIB_{params.UIB:.2f}_M_{grid.M}_N_{params.N}.hdf5'

        with h5py.File(filename, 'r') as f:
            dJU_values = f['dJU_values'][:]
            mu_qcorr = f['mu_qcorr'][:]
        return dJU_values, mu_qcorr

    def read_from_hdf5_perturbative(self, UIB, muU, M, N):
        filename = f'../data/Perturbative_UIB_{UIB:.2f}_Mu_{muU}_M_{M}_N_{N}.hdf5'
        with h5py.File(filename, 'r') as f:
            dJU_values = f['dJU_values'][:]
            # en_values = f['en_values'][:]
            omega0s = f['omega0s'][:]
            omega1s = f['omega1s'][:]
            omega2s = f['omega2s'][:]
        return dJU_values, omega0s, omega1s, omega2s

    def read_from_hdf5_all_fixed_chemical_potential(self, UIB, muU, M, N):
        filename = f'../data/spectral_funcs_UIB_{UIB:.2f}_Mu_{muU}_M_{M}_N_{N}.hdf5'
        with h5py.File(filename, 'r') as f:
            dJU_values = f['dJU_values'][:]
            # en_values = f['en_values'][:]
            omega0s = f['omega0s'][:]
            omega1s = f['omega1s'][:]
            omega2s = f['omega2s'][:]

            T11_values = []
            T12_values = []
            T21_values = []
            T22_values = []
            T22_SE_values = []
            SE_SI = []

            for dJU in dJU_values:
                group = f[f'dJU_{dJU}']
                T11_values.append(group['T11'][:])
                T12_values.append(group['T12'][:])
                T21_values.append(group['T21'][:])
                T22_values.append(group['T22'][:])
                T22_SE_values.append(group['T22_SE'][:])
                SE_SI.append(group['SE_SI'][:])
                en_values = group['en'][:]

        return dJU_values, en_values, omega0s, omega1s, omega2s, T11_values, T12_values, T21_values, T22_values, T22_SE_values, SE_SI
    
    def read_from_hdf5_all_fixed_density(self, UIB, n0, M, N):
        filename = f'../data/spectral_funcs_UIB_{UIB:.2f}_n0_{n0}_M_{M}_N_{N}.hdf5'
        with h5py.File(filename, 'r') as f:
            dJU_values = f['dJU_values'][:]
            # en_values = f['en_values'][:]
            omega0s = f['omega0s'][:]
            omega1s = f['omega1s'][:]
            omega2s = f['omega2s'][:]

            T11_values = []
            T12_values = []
            T21_values = []
            T22_values = []
            T22_SE_values = []
            SE_SI = []

            for dJU in dJU_values:
                group = f[f'dJU_{dJU}']
                T11_values.append(group['T11'][:])
                T12_values.append(group['T12'][:])
                T21_values.append(group['T21'][:])
                T22_values.append(group['T22'][:])
                T22_SE_values.append(group['T22_SE'][:])
                SE_SI.append(group['SE_SI'][:])
                en_values = group['en'][:]

        return dJU_values, en_values, omega0s, omega1s, omega2s, T11_values, T12_values, T21_values, T22_values, T22_SE_values, SE_SI
    # def read_from_hdf5(self, UIB, muU, M, N):
    #     filename = f'../data/spectral_funcs_UIB_{UIB:.2f}_Mu_{muU}_M_{M}_N_{N}.hdf5'

    #     with h5py.File(filename, 'r') as f:
    #         # Read omega0s, omega1s, omega2s
    #         omega0s = f['omega0s'][:]
    #         omega1s = f['omega1s'][:]
    #         omega2s = f['omega2s'][:]
    #         dJUs = f['dJU_values'][:]
    #         # en_vector = f['en_values'][:]
    #         # T11_read = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
    #         # T12_read = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
    #         # T21_read = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
    #         # T22_read = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
    #         # T22_SE_read = np.zeros((len(dJUs), len(en_vector)), dtype=np.complex128)
    #         # Read the data for each dJU value
    #         dJU_values = [key for key in f.keys() if key.startswith('dJU_')]
    #         data = {}
    #         for dJU in dJU_values:
    #             counter = 0
    #             group = f'dJU_{dJU}'
    #             data[dJU] = {
    #                 'T11': group['T11'][:],
    #                 'T12': group['T12'][:],
    #                 'T21': group['T21'][:],
    #                 'T22': group['T22'][:],
    #                 'T22_SE': group['T22_SE'][:],
    #                 'SE_SI': group['SE_SI'][:],
    #                 'en': group['en'][:]
    #             }
    #             # T11_read[counter, :] = group['T11'][:]
    #             # T12_read[counter, :] = group['T12'][:]
    #             # T21_read[counter, :] = group['T21'][:]
    #             # T22_read[counter, :] = group['T22'][:]
    #             # T22_SE_read[counter, :] = group['T22_SE'][:]
    #             counter += 1
    #             # en_vector = group['en'][:]
    #         # Read the attributes
    #         UIB = f.attrs['UIB']
    #         Lx = f.attrs['Lx']
    #         Ly = f.attrs['Ly']
    #         muU = f.attrs['muU']

    #     return omega0s, omega1s, omega2s, data, UIB, Lx, Ly, muU
    
