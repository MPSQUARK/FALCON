from astropy.io import fits
import h5py
import numpy as np

target_dtype = '<f4'
comp_type = 'gzip'
comp_level = 9

ngp = fits.getdata('DustMaps/SFD_dust_4096_ngp.fits')
sgp = fits.getdata('DustMaps/SFD_dust_4096_sgp.fits')

# Convert to little endian
ngp = ngp.astype(target_dtype, copy=False)
sgp = sgp.astype(target_dtype, copy=False)

with h5py.File('DustMaps/dust.h5', 'w') as f:
	f.create_dataset('ngp', data=ngp, compression=comp_type, chunks=True, compression_opts=comp_level, dtype=target_dtype)
	f.create_dataset('sgp', data=sgp, compression=comp_type, chunks=True, compression_opts=comp_level, dtype=target_dtype)