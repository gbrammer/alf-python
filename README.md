# alf-python
Python bindings to Charlie Conroy's Alf code: https://github.com/cconroy20/alf.

## Installation

(Install open-mpi with Homebrew.)

```bash
python setup.py install
```

## Demo

```python
import matplotlib.pyplot as plt
import numpy as np

from alf.alf import Alf

# Initialize the Fortran functions
sps = Alf()
# Alf: CALL SETUP()

sp = sps.get_model(in_place=False, logage=0.1)
plt.plot(sps.wave, sp)

sp = sps.get_model(in_place=False, logage=0.5, sigma=150)
plt.plot(sps.wave, sp)

sps.show_params()
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

# 0: full fit, 1: limited abundances, 2: simple fit
print(sps.fit_type) 
sps.set_fit_type(1)

# Show variation with, e.g., mgh
sps.initialize_defaults()
m = sps.get_model(in_place=False)

for mgh in np.arange(-0.8,0.11,0.1):
    sps.set_param(mgh=mgh)
    m_i = sps.get_model(in_place=False)
    plt.plot(sps.wave, m/m_i, alpha=0.5, label='{0:.1f}'.format(mgh))
    
```
