# alf-python
Python bindings to Charlie Conroy's Alf code: https://github.com/cconroy20/alf.

## Installation

Hacky way to deal with f2py and Alf compiled with MPIFORT:

```
$ python setup.py install
$ sh force_mpifort.sh
```

The latter script breaks ```python setup.py install``` intentionally and then runs the last compliation command with mpifort rather than gfortran.  Finally the script copies the newly-created `so` file to the PATH of the `alf` module (hardcoded for an OSX / Python 3.5 anaconda distro).

## Demo

```python
import matplotlib.pyplot as plt
import alf.alf

alf_sps = alf.alf.Alf()
# Alf: CALL SETUP()

sp = alf_sps.get_model(in_place=False, logage=0.1)
plt.plot(alf_sps.wave, sp)

sp = alf_sps.get_model(in_place=False, logage=0.5, sigma=150)
plt.plot(alf_sps.wave, sp)

alf_sps.show_params()
           velz =     0.00
          sigma =   150.00
         logage =     0.50
             zh =     0.00
            feh =     0.00
             ah =     0.00
             ch =     0.00
             nh =     0.00
            nah =     0.00
            mgh =     0.00
            sih =     0.00
             kh =     0.00
            cah =     0.00
            tih =     0.00
             vh =     0.00
            crh =     0.00
            mnh =     0.00
            coh =     0.00
            nih =     0.00
            cuh =     0.00
            srh =     0.00
            bah =     0.00
            euh =     0.00
           teff =     0.00
           imf1 =     1.30
           imf2 =     2.30
          logfy =    -4.00
         sigma2 =     0.00
          velz2 =     0.00
         logm7g =    -4.00
        hotteff =    20.00
         loghot =    -4.00
      fy_logage =     0.30
       logtrans =    -4.00
    logemline_h =    -4.00
 logemline_oiii =    -4.00
  logemline_sii =    -4.00
   logemline_ni =    -4.00
  logemline_nii =    -4.00
         jitter =     1.00
           imf3 =     2.00
         logsky =    -4.00
           imf4 =     0.00
             h3 =     0.00
             h4 =     0.00

```
