import numpy as np
from ._alf import driver

ALF_PARAMS = [p.strip() for p in """
velz          
sigma         
logage        
zh            
feh           
ah            
ch            
nh            
nah           
mgh           
sih           
kh            
cah           
tih           
vh            
crh           
mnh           
coh           
nih           
cuh           
srh           
bah           
euh           
teff          
imf1          
imf2          
logfy         
sigma2        
velz2         
logm7g        
hotteff       
loghot        
fy_logage     
logtrans      
logemline_h   
logemline_oiii
logemline_sii 
logemline_ni  
logemline_nii 
jitter        
imf3          
logsky        
imf4          
h3            
h4""".split()]

pidx = {}
for i, k in enumerate(ALF_PARAMS):
    pidx[k] = i

class Alf(object):
    def __init__(self, verbose=True, ldeg=400, mask_lines=False):
                
        if not bool(driver.is_setup):
            if verbose:
                print('Alf: CALL SETUP()')
        
            driver.setup_alf()
        
        self.nl = driver.get_nspec()
        self.npar = driver.get_npar() 
        self.nfil = driver.get_nfil() 
        
        # Polynomial order for every `ldeg` angstroms of available spectrum
        self.ldeg = ldeg
        
        self.MASK_LINES = mask_lines
        
        self.default_params = driver.get_default_parameters(self.npar)
        self.wave = driver.get_grid_lam(self.nl)
        self.params = self.default_params*1
        self.get_model()
        
    def get_model(self, in_place=True, **kwargs):
        self.set_param(**kwargs)
        # These must be in different order than in alf.f90, don't know why
        sp = driver.get_spec(self.nl, self.params, self.npar)
        if in_place:
            self.spec = sp
        else:
            return sp
    
    def get_M2L(self, in_place=True, **kwargs):
        """Get line-free model and M/L ratios"""
        self.set_param(**kwargs)
        # These must be in different order than in alf.f90, don't know why
        m2l = driver.get_m2l(self.nl, self.nfil, self.params, self.npar)

        return m2l   
                 
    def set_param(self, **kwargs):
        for k in kwargs:
            if k not in ALF_PARAMS:
                print('Paramter "{0}" not valid.'.format(k))
            else:
                self.params[pidx[k]] = kwargs[k] 
    
    def info(self, verbose=True):
        lines = []
        for k in ALF_PARAMS:
            lines.append('{0:>15s} = {1:8.2f}'.format(k, self.params[pidx[k]]))
        
        if verbose:
            print('\n'.join(lines))
        else:
            return lines
    
    @staticmethod
    def get_default_priors(param_names):
        from scipy.stats import uniform
        from collections import OrderedDict
        
        prior_limits = {'velz':(-1000,1000),          
                         'sigma':(0,1000),       
                         'logage':(0.1, 1.0),    
                         'zh':  (-1, 0.2),          
                         'feh': (-0.25, 0.25),          
                         'ah':  (-0.25, 0.25),          
                         'ch':  (-0.25, 0.25),          
                         'nh':  (-0.25, 0.25),          
                         'nah': (-0.25, 0.25),          
                         'mgh': (-0.25, 0.25),          
                         'sih': (-0.25, 0.25),          
                         'kh':  (-0.25, 0.25),          
                         'cah': (-0.25, 0.25),          
                         'tih': (-0.25, 0.25),          
                         'vh':  (-0.25, 0.25),          
                         'crh': (-0.25, 0.25),          
                         'mnh': (-0.25, 0.25),          
                         'coh': (-0.25, 0.25),          
                         'nih': (-0.25, 0.25),          
                         'cuh': (-0.25, 0.25),          
                         'srh': (-0.25, 0.25),          
                         'bah': (-0.25, 0.25),          
                         'euh': (-0.25, 0.25),
                         'logemline_h':   (-4,-2),
                         'logemline_oiii':(-4,-2),
                         'logemline_sii': (-4,-2),
                         'logemline_ni':  (-4,-2),
                         'logemline_nii': (-4,-2),
                         'slope':(-20,20)}
        
        prior = OrderedDict()
        for p in param_names:
            if p in prior_limits:
                lim = prior_limits[p]
                prior[p] = uniform(loc=lim[0], scale=lim[1]-lim[0])
            else:
                prior[p] = uniform(loc=-1.e10, scale=2.e10)
        
        return prior
    
    @staticmethod
    def evaluate_prior(param_values, param_names, prior):
        lnp_prior = 0
        for ip, p in enumerate(param_names):
            if p in prior:
                lnp_prior += prior[p].logpdf(param_values[ip])
        
        return lnp_prior
        
    @staticmethod
    def get_parameter_scaling(param_names):
        default_scale = {'velz':5,          
                         'sigma':50,       
                         'logage':0.1,    
                         'zh': 0.1,          
                         'feh':0.1,          
                         'ah': 0.1,          
                         'ch': 0.1,          
                         'nh': 0.1,          
                         'nah':0.1,          
                         'mgh':0.1,          
                         'sih':0.1,          
                         'kh': 0.1,          
                         'cah':0.1,          
                         'tih':0.1,          
                         'vh': 0.1,          
                         'crh':0.1,          
                         'mnh':0.1,          
                         'coh':0.1,          
                         'nih':0.1,          
                         'cuh':0.1,          
                         'srh':0.1,          
                         'bah':0.1,          
                         'euh':0.1,
                         'logemline_h':2,
                         'logemline_oiii':2,
                         'logemline_sii':2,
                         'logemline_ni':2,
                         'logemline_nii':2}
        
        scale = np.ones(len(param_names))
        for ip, p in enumerate(param_names):
            if p in default_scale:
                scale[ip] = default_scale[p]
                
        #scale = np.array([default_scale[p] for p in param_names])
        return scale
        
    @staticmethod
    def _obj_fit_alf(fit_params, alf_data, sps, param_names, scale, retval):
        """
        Objective function for fitting Alf spectra
        
        alf_data : np.ndarray([rest_wave, flux, err, mask, inst])
        
        sps : Alf instance
        
        param_names : Alf paramter names corresponding to `fit_params`
        
        scale : parameter scaling factors, array or float
            Rescaling helps scipy.optimize.minimize functions
        
        retval : control what is returned
            'model': mask, model, poly arrays
            'resid': resdiuals / err for scipy.optimize.least_squares
            'lnp':   log probability for EMCEE
            'chi2':  chi-squared
            
        """
        #velz, sigma, logage, cah = fit_params
        if isinstance(scale, (list, tuple, np.ndarray)):
            scale_array = scale
        else:
            scale_array = np.ones(len(fit_params))

        param_keys = {}
        for i, k in enumerate(param_names):
            param_keys[k] = fit_params[i]*scale_array[i]

        if 'sigma' in param_names:  
            if param_keys['sigma'] > 9999:
                if retval == 'lnp':
                    return np.inf
                elif retval == 'resid':
                    return np.ones(alf_data.shape[1])*np.inf
                else:
                    return -np.inf

        model = sps.get_model(in_place=False, **param_keys)

        if 'velz' in param_names:
            oneplusz = (1-param_keys['velz']/CLIGHT*1E5)
        else:
            oneplusz = 1.

        mask = (alf_data[0,:] > 3700) & (alf_data[0,:] < 8500)
        if sps.MASK_LINES:
            for linew in [3727., 4102., 4341., 4862., 4960., 5008., 5203., 6549., 6564., 6585., 6718., 6732.]:
                line_mask = np.abs(alf_data[0,:]-linew)/linew*CLIGHT/1.e5 > 500.
                mask &= line_mask

        modelz = np.interp(alf_data[0,mask], sps.wave/oneplusz, model, left=0, right=0)

        deg = int((alf_data[0,mask].max()-alf_data[0,mask].min())/sps.ldeg)
        poly_coeff = np.polyfit(alf_data[0,mask], alf_data[1,mask]/modelz, deg, w=1./alf_data[2,mask])

        if retval == 'model':
            modelz = np.interp(alf_data[0,:], sps.wave/oneplusz, model, left=0, right=0)
            scale_poly = np.polyval(poly_coeff, alf_data[0,:])
            return mask, modelz, scale_poly

        resid = (alf_data[1,mask]-modelz*np.polyval(poly_coeff, alf_data[0,mask]))/alf_data[2,mask]

        if retval == 'resid':
            return resid

        chi2 = (resid**2).sum()

        if retval == 'lnp':        
            print(fit_params, chi2)
            if not np.isfinite(chi2):
                return -np.inf
            else:
                return -0.5*chi2

        #print(fit_params, chi2)

        return chi2

    @classmethod
    def array_to_dict(self, params, param_names=[]):
        from collections import OrderedDict
        p_dict = OrderedDict()
        for i, p in enumerate(param_names):
            p_dict[p] = params[i]
        
        return p_dict

def show_jac(res, wave, param_names):
    """
    Make a plot of the Jacobian result from a least_squares 'lm' 
    fit to show the parameter derivatives as a function of wavelength
    """
    import matplotlib.pyplot as plt
    N = len(param_names)
    fig = plt.figure(figsize=[8,1*N])
    for i in range(N):
        ax = fig.add_subplot(N,1,i+1)
        ax.plot(wave, res.jac[:,i], color='k')
        if i < (N-1):
            ax.set_xticklabels([])
        
        ax.grid()
        ax.set_ylabel(param_names[i])
        ax.set_ylim(-9,9)
    
    fig.tight_layout(pad=0.1)
    
    