
![Logo](https://suit.iucaa.in/sites/default/files/top_banner_compressed_2_1.png)


# PSF Generator from Zernike Coefficients

Generate PSFs from Zernike Coefficients with `psf_generator_cartesian.py`. Also, generate wavefront map for unit circle.

Steps:
- Wavefront map is generated from a list of Zernike Coefficients.
- We use [Piston, TiltX, TiltY, Defocus, AstigX, AstigY, PrimarySph, TrefoilX, TrefoilY]
- PSF is given by the square of the absolute value of the FFT of complex pupil function.
- The complex pupil function is generated with the Zernike Coefficients.
## Acknowledgements

 - [IUCAA, Pune](https://www.iucaa.in)
 - [ISRO, Aditya-L1](https://www.isro.gov.in/Aditya_L1.html)

## Authors

- [@janmejoysarkar](https://github.com/janmejoysarkar)
