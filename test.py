from pygsfit_cp import *
gs = pygsfit_cp(filename='/Users/walterwei/Downloads/20220511/slf_final_XX_t19_allbd.fits', out_dir='/Users/walterwei/Downloads/20220511/gsfit_test')
print('The EOVSA image cube is loaded, which is in size of {0}(freqs) * {1}pix * {2}pix'.format(*gs.flux_data.shape))
gs.ninput = np.array([7, 0, 1, 30, 1, 1], dtype='int32') #Nparms, Angular_mode,Npix, Nfreq(will be replaced) fitting mode, stokes
gs.rinput = np.array([0.17, 1e-6, 1.0, 4.0, 8.0, 0.015], dtype='float64')  # real_input
gs.initial_parguess = np.array([ # *Value, *lower boundary * upper boundary
                [10.0, 0.0001, 2000.0], #n_nth;    1d7 cm^-3
                [4.0, 0.01, 30.0], #B;    1d2G
                [60.0, 22.0, 87.0], #theta;    deg
                [10.0, 0.01, 600.0], #n_th;    1d9 cm^-3
                [4.5, 1.6, 10.0], #Delta;    No
                [5.0, 0.1, 10.0], #E_max;    MeV
                [5.0, 1.5, 60.0], #T_e;    MK
            ], dtype=np.float64)
gs.update_flux_threshold_mask(threshold=8)
gs.do_fit(mode='batch')