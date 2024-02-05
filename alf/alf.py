from collections import OrderedDict

import numpy as np
from ._alf import driver
from astropy.cosmology import Planck15

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
logemline_h   
logemline_oii
logemline_oiii
logemline_sii 
logemline_ni  
logemline_nii 
logtrans      
jitter        
logsky        
imf3          
imf4          
h3            
h4
""".split()]

EMISSION_LINE_PARAMS = ['logemline_h', 'logemline_oii', 
                        'logemline_oiii', 'logemline_sii', 
                        'logemline_ni', 'logemline_nii']

SECOND_COMPONENT_PARAMS = ['logfy', 'fy_logage']
SECOND_COMPONENT_VEL_PARAMS = ['sigma2','velz2']
                 
PRIOR_LIMITS = {'velz':[-500,500],
                 'sigma':[10,500],
                 'velz_line':[-500,500],
                 'sigma_line':[10,500],
                 'velz2':[-1000,1000],
                 'sigma2':[10,1000],
                 'logage':[np.log10(0.5), np.log10(14)], 
                 'zh':  [-2.0, 0.25],
                 'feh': [-0.3, 0.5],
                 'ah':  [-0.3, 0.5],
                 'ch':  [-0.3, 0.5],
                 'nh':  [-0.3, 1.0],
                 'nah': [-0.3, 1.0],
                 'mgh': [-0.3, 0.5],
                 'sih': [-0.3, 0.5],
                 'kh':  [-0.01, 0.01], # used as a proxy for grism shift since doesn't (normally) affect spectrum          
                 'cah': [-0.3, 0.5],
                 'tih': [-0.3, 0.5],
                 'vh':  [-0.3, 0.5],
                 'crh': [-0.3, 0.5],
                 'mnh': [-0.3, 0.5],
                 'coh': [-0.3, 0.5],
                 'nih': [-0.3, 0.5],
                 'cuh': [-0.3, 0.5],
                 'srh': [-0.3, 0.5],
                 'bah': [-0.6, 0.5],
                 'euh': [-0.6, 0.5],
                 'teff': [-50, 50],
                 'logfy': [-6, -0.05],
                 'fy_logage': [np.log10(0.5), -0.15],
                 'logm7g': [-6.0, -1.0],
                 'hotteff': [8.0, 30.],
                 'loghot': [-6., -1],
                 'logtrans': [-6.0, 1.],
                 'jitter': [0.1, 10.],
                 'logemline_h':   [-6.,2],
                 'logemline_oii':[-6,2],
                 'logemline_oiii':[-6,2],
                 'logemline_sii': [-6,2],
                 'logemline_ni':  [-6,2],
                 'logemline_nii': [-6,2],
                 'logsky': [-9,3],
                 'slope':[-20,20],
                 'dzgrism':[-0.01, 0.01],
                 'h3':[-0.4, 0.4],
                 'h4':[-0.4, 0.4]}


pidx = {}
for i, k in enumerate(ALF_PARAMS):
    pidx[k] = i


class Alf(object):
    def __init__(self, fit_type=1, fit_two_ages=0, maskem=1, fit_hermite=0, velbroad_simple=0, verbose=True, ldeg=400, mwimf=1, imf_type=1, mask_lines=False):
        """Container for the ``alf`` driver
        """
        if not bool(driver.is_setup):
            if verbose:
                print('Alf: CALL SETUP()')
        
            driver.setup_alf()
        
        self.nl = driver.get_nspec()
        self.npar = driver.get_npar() 
        self.nfil = driver.get_nfil() 
        self.nage = driver.get_nage() 
        self.ssp_logage = driver.get_ssp_logagegrid(self.nage)
        
        self.set_fit_type(fit_type)
        self.set_fit_two_ages(fit_two_ages)
        
        self.set_imf(mwimf=mwimf, imf_type=imf_type)
        self.set_maskem(maskem)
        self.set_fit_hermite(fit_hermite)
        
        # Polynomial order for every `ldeg` angstroms of available spectrum
        self.ldeg = ldeg
        
        self.MASK_LINES = mask_lines
        
        self.default_params = driver.get_default_parameters(self.npar)
        
        self.wave = driver.get_grid_lam(self.nl)
        self.set_wave_limits(self.wave[0], self.wave[-1])
        
        self.params = self.default_params*1
        self.param_names = ALF_PARAMS
        self.fit_params = [p for p in self.param_names]
        self.get_model()
        
        self.dummy_params = {}
        
        self.prior = self.get_default_priors(list(PRIOR_LIMITS.keys()), 
                                            limits=PRIOR_LIMITS)


    def set_fit_type(self, fit_type):
        """
        Set fit type: 
            !0: fit the full model (IMF, all abundances, nuisance params, etc)
            !1: only fit velz, sigma, SSP age, Z, Fe,C,N,O,Mg,Si,Ca,Ti,Na
            !2: only fit velz, sigma, SSP age, Z
        """
        driver.set_fit_type(int(fit_type))
        self.fit_type = fit_type


    def set_wave_limits(self, il1, il2):
        """
        Set wavelength limits, which will only be used for limiting the 
        wavelengths over which VELBROAD is run when `self.fit_hermite==1`.
        """
        driver.set_wave_limits(il1, il2)
        self.l1 = il1
        self.l2 = il2


    def set_maskem(self, maskem):
        """
        Set maskem for emission lines: 
            !0: Include emission lines in model
        """
        driver.set_maskem(int(maskem))
        self.maskem = maskem


    def set_fit_hermite(self, fit_hermite):
        """
        Set fit_hermite and velbroad_simple for using hermite coefficients
        """
        driver.set_fit_hermite(int(fit_hermite))
        driver.set_velbroad_simple(int(fit_hermite))
        self.velbroad_simple = fit_hermite
        self.fit_hermite = fit_hermite


    def set_fit_two_ages(self, fit_two_ages, age_grid=[0.5, 1, 14]):
        """
        Set fit_two_ages flag.
        
        age_split: [tmin, tsplit, tmax], None
            If not None, reset priors on 'logage' and 'fy_logage' such that 
            they don't overlap:
            
                fy_logage = uniform[min..split]
                logage = uniform[split..max]
                
               
        """
        from scipy.stats import uniform
        driver.set_fit_two_ages(int(fit_two_ages))
        
        if fit_two_ages > 0:
            self.set_fit_type(0)
            
            if age_grid is not None:
                tmin, tsplit, tmax = np.log10(age_grid)
                self.prior['fy_logage'] = uniform(loc=tmin, scale=tsplit-tmin)
                self.prior['logage'] = uniform(loc=tsplit, scale=tmax-tsplit)
                    
        self.fit_two_ages = fit_two_ages
    
    def set_imf(self, mwimf=1, imf_type=1):
        """
        Set IMF params
        """
        driver.set_imf(int(mwimf), int(imf_type))
        self.mwimf = mwimf
        self.imf_type = imf_type


    def get_model(self, in_place=True, **kwargs):
        """Generate a model spectrum
        
        Parameters
        ----------
        in_place : bool
            Put result in ``spec`` attribute
        
        kwargs : dict
            Parameters to modify at runtime, passed to ``set_param``
        
        """
        self.set_param(**kwargs)
        # These must be in different order than in alf.f90, don't know why
        if self.mwimf:
            sp = driver.get_spec_mw(self.nl, self.params, self.npar)
        else:
            sp = driver.get_spec(self.nl, self.params, self.npar)
            
        if in_place:
            self.spec = sp
        else:
            return sp


    def interp_model(self, z=0, wave=None, interp_func=np.interp, **kwargs):
        """
        Interpolate Alf model on a new wavelength grid
        """
        self.set_param(**kwargs)
        sp = driver.get_spec(self.nl, self.params, self.npar)

        ## Don't do any interpolation
        if wave is None:
            return sp
        else:
            interp_sp = interp_func(wave, self.wave*(1+z), sp)
            return interp_sp


    def spline_matrix(self, z=0, wave=None, knots=None, interp_func=np.interp, **kwargs):
        """
        Get spline terms
        """
        from patsy import dmatrix
        
        interp_sp = self.interp_model(z=z, wave=wave, interp_func=interp_func, **kwargs)
        if wave is None:
            wave = self.wave
          
        if knots is None:
            knots = np.percentile(wave, [20,40,60,80])
        
        data = {'w':wave, 'k':knots, 'model':interp_sp}
        _Xsp = dmatrix('model*bs(w, df=None, knots=k, degree=3) - 1', data)
        return _Xsp


    def get_M2L(self, in_place=True, **kwargs):
        """Get line-free model and M/L ratios"""
        self.set_param(**kwargs)
        # These must be in different order than in alf.f90, don't know why
        m2l = driver.get_m2l(self.nl, self.nfil, self.params, self.npar)

        return m2l   


    def set_param(self, verbose=False, **kwargs):
        """Update parameters in ``params`` dictionary
        """
        for k in kwargs:
            if k in ALF_PARAMS:
                self.params[pidx[k]] = kwargs[k] 
            elif k in self.dummy_params:
                self.dummy_params[k] = kwargs[k]
            else:
                print('Parameter "{0}" not valid.'.format(k))
    
    @property
    def param_dict(self):
        """
        Dictionary of available parameters for a given `fit_type`.
        """
        from collections import OrderedDict
        d = OrderedDict()
        
        # Number of properties for a given fit_type
        Np = [46,14,4][self.fit_type]
        allowed_params = self.param_names[:Np]
        
        if self.fit_hermite > 0:
            allowed_params += ['h3','h4']
        
        if self.fit_two_ages > 0:
            allowed_params += ['logfy', 'fy_logage']
            
        if (self.maskem == 0) & (self.fit_type == 0):
            allowed_params += EMISSION_LINE_PARAMS
                                       
        for i, (name, p) in enumerate(zip(self.param_names, self.params)):            
            if (name in self.fit_params) & (name in allowed_params):
                d[name] = p
        
        for k in self.dummy_params:
            d[k] = self.dummy_params[k]
            
        return d


    def __setitem__(self, p, value):
        self.set_param(**{p:value})


    def __getitem__(self, p):
        ix = self.param_names.index(p)
        return self.params[ix]

        
    def info(self, verbose=True):
        """Print parameter information
        """
        lines = []
        for k in self.param_names:
            lines.append('{0:>15s} = {1:8.2f}'.format(k, self.params[pidx[k]]))
        
        if verbose:
            print('\n'.join(lines))
        else:
            return lines


    def initialize_defaults(self):
        """Set default parameter values
        """
        self.set_param(**self.array_to_dict(self.default_params, ALF_PARAMS))

    
    @classmethod
    def array_to_dict(self, params, param_names=[]):
        """Create parameter dictionary from array list
        """
        from collections import OrderedDict
        p_dict = OrderedDict()
        ndim = params.ndim
        for i, p in enumerate(param_names):
            if ndim == 1:
                p_dict[p] = params[i]
            elif ndim == 2:
                p_dict[p] = params[:,i]
            elif ndim == 3:
                p_dict[p] = params[:,:,i]
            else:
                p_dict[p] = None
                
        return p_dict


    def get_defaults(self, param_names, **kwargs):
        """
        Get default parameter values
        """
        pdict = self.array_to_dict(self.default_params, ALF_PARAMS)
        plist = [pdict[p] for p in param_names]
        for k in kwargs:
            if k in param_names:
                plist[param_names.index(k)] = kwargs[k]
                
        return np.array(plist)


    @staticmethod
    def get_default_priors(param_names, limits=PRIOR_LIMITS):
        """
        Get default prior ranges
        """
        from scipy.stats import uniform
        from collections import OrderedDict
                
        prior = OrderedDict()
        for p in param_names:
            if p in limits:
                lim = limits[p]
                prior[p] = uniform(loc=lim[0], scale=lim[1]-lim[0])
            else:
                prior[p] = uniform(loc=-1.e10, scale=2.e10)
        
        return prior


    def evaluate_prior(self):
        """
        Evaluate prior based on param_dict
        """
        pd = self.param_dict
        lnp = 0
        for p in pd:
            lnp += self.prior[p].logpdf(pd[p])
        
        return lnp


    def get_prior_limits(self, param_names):
        """Get prior ranges for a list of parameter names
        """
        limits = [PRIOR_LIMITS[p] for p in param_names]
        return limits


    @staticmethod
    def get_parameter_scaling(param_names):
        """
        Rescaling  values for a list of parameters
        """
        default_scale = {'velz':50,
                         'sigma':50,
                         'velz_line':50,
                         'sigma_line':50,
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
                         'logemline_h':   0.4,
                         'logemline_oii': 0.4,
                         'logemline_oiii':0.4,
                         'logemline_sii': 0.4,
                         'logemline_ni':  0.4,
                         'logemline_nii': 0.4,
                         'dzgrism':0.001}
        
        scale = np.ones(len(param_names))
        for ip, p in enumerate(param_names):
            if p in default_scale:
                scale[ip] = default_scale[p]
                
        #scale = np.array([default_scale[p] for p in param_names])
        return scale


    def get_mass_weighted_age(self):
        """Calculate mass-weighted age of the current stellar population
        """
        p0 = self.array_to_dict(self.params, self.param_names)
        fy = 10**(p0['logfy'])
        mass_age = fy*10**p0['fy_logage']+(1-fy)*10**(p0['logage'])
        return mass_age


    def abundance_correct(self, pdict=None, s07=False, b14=False, m11=True):
        """
        Abundance corrections
        """
        from scipy import constants, interpolate

        p0 = self.array_to_dict(self.params, self.param_names)
        
        if pdict is None:
            pdict = self.array_to_dict(self.params, self.param_names)
        
        # Find what parameters are specified and what are part of the `self`
        # model.
        for c in self.param_names:
            if c.endswith('h') & (~c.startswith('logemline')):
                if c not in pdict:
                    #print(c)
                    pdict[c] = pdict['logage']*0+p0[c]
                            
        ## From read_alf.py
        # Correction factros from Schiavon 2007, Table 6
        # NOTE: Forcing factors to be 0 for [Fe/H]=0.0,0.2
        lib_feh = [-1.6, -1.4, -1.2, -1.0, -0.8,
                   -0.6, -0.4, -0.2, 0.0, 0.2]
        lib_ofe = [0.6, 0.5, 0.5, 0.4, 0.3, 0.2,
                   0.2, 0.1, 0.0, 0.0]

        if s07:
            #Schiavon 2007
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.29,
                        0.20, 0.13, 0.08, 0.05, 0.04]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.20,
                        0.12, 0.06, 0.02, 0.0, 0.0]
        elif b14:
            # Fitted from Bensby+ 2014
            lib_mgfe = [0.4 , 0.4, 0.4, 0.38, 0.37,
                        0.27, 0.21, 0.12, 0.05, 0.0]
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]
        elif m11 or (b14 is False and s07 is False):
            # Fitted to Milone+ 2011 HR MILES stars
            lib_mgfe = [0.4, 0.4, 0.4, 0.4, 0.34, 0.22,
                        0.14, 0.11, 0.05, 0.04]
            # from B14
            lib_cafe = [0.32, 0.3, 0.28, 0.26, 0.26,
                        0.17, 0.12, 0.06, 0.0, 0.0]

        # In ALF the oxygen abundance is used
        # a proxy for alpha abundance
        del_alfe = interpolate.interp1d(lib_feh, lib_ofe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_mgfe = interpolate.interp1d(lib_feh, lib_mgfe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')
        del_cafe = interpolate.interp1d(lib_feh, lib_cafe,
                                        kind='linear',
                                        bounds_error=False,
                                        fill_value='extrapolate')

        #zh = np.where(self.labels == 'zH')
        if np.isscalar(pdict['zh']):
            al_corr = float(del_alfe(pdict['zh']))
            mg_corr = float(del_mgfe(pdict['zh']))
            ca_corr = float(del_cafe(pdict['zh']))
        else:
            al_corr = del_alfe(pdict['zh'])
            mg_corr = del_mgfe(pdict['zh'])
            ca_corr = del_cafe(pdict['zh'])
            
        group0 = ['a', 'Mg']
        
        # Assuming Ca~Ti~Si
        group1 = ['Ca', 'Ti', 'Si']

        # These elements seem to show no net enhancemnt
        # at low metallicity
        group2 = ['C', 'Ca', 'N', 'Cr', 'Ni', 'Na']

        # These elements we haven't yet quantified
        group3 = ['Ba', 'Eu', 'Sr', 'Cu', 'Co', 'K', 'V', 'Mn']
        
        # alpha
        pdict['aFe'] = pdict['ah'] - pdict['feh'] + al_corr
        
        # Mg
        pdict['MgFe'] = pdict['mgh'] - pdict['feh'] + mg_corr
        
        # Ca-like
        for elem in group1:
            ci = elem.lower()+'h'
            pdict['{0}Fe'.format(elem)] = pdict[ci] - pdict['feh'] + ca_corr
        
        # Others
        for elem in np.hstack([group2, group3]):
            ci = elem.lower()+'h'
            pdict['{0}Fe'.format(elem)] = pdict[ci] - pdict['feh']
        
        return pdict


    def get_stellar_mass(self, z=0, rnorm=1., cosmo=Planck15):
        """
        Get M/L, L, stellar mass, as measured in R-band
        
        Returns: 
            M/Lr, log10(Lr), log10(stellar_mass)
            
        """
        from sedpy.observate import Filter
        import astropy.units as u
        import astropy.constants as const
        
        MLr, MLi, MLk = self.get_M2L()
        #plt.plot(self.wave, self.spec); print(MLr)
        
        rfilt = Filter('sdss_r0')
        
        norm_spec = self.get_model(in_place=False)*rnorm
        
        #rband_Jy = rfilt.obj_counts(self.wave, norm_spec)/rfilt.ab_zero_counts*3631
        #rband_flam = rband_Jy/(3.34e4*rfilt.wave_pivot**2)#*u.erg/u.second/u.cm**2/u.AA
        #dL = Planck15.luminosity_distance(z)
        
        Mr = rfilt.ab_mag(self.wave, norm_spec) - cosmo.distmod(z).value
        Lr = 10**(-0.4*(Mr-rfilt.solar_ab_mag))
        #Lr = (rband_flam*4*np.pi*dL.to(u.cm).value**2)/3.828e+33*rfilt.wave_pivot*(1+z)
        stellar_mass = Lr*MLr
        return MLr, np.log10(Lr), np.log10(stellar_mass)


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
            oneplusz = (1 - param_keys['velz']/CLIGHT*1E5)
        else:
            oneplusz = 1.

        mask = (alf_data[0,:] > 3700) & (alf_data[0,:] < 8500)
        if sps.MASK_LINES:
            for linew in [3727., 4102., 4341., 4862., 4960., 5008.,
                          5203., 6549., 6564., 6585., 6718., 6732.]:
                line_mask = np.abs(alf_data[0,:]-linew)/linew*CLIGHT/1.e5 > 500.
                mask &= line_mask

        modelz = np.interp(alf_data[0,mask], sps.wave/oneplusz, model, left=0, right=0)

        deg = int((alf_data[0,mask].max()-alf_data[0,mask].min())/sps.ldeg)
        poly_coeff = np.polyfit(alf_data[0,mask],
                                alf_data[1,mask]/modelz, deg,
                                w=1./alf_data[2,mask])

        if retval == 'model':
            modelz = np.interp(alf_data[0,:], sps.wave/oneplusz, model, left=0, right=0)
            scale_poly = np.polyval(poly_coeff, alf_data[0,:])
            return mask, modelz, scale_poly

        resid = (alf_data[1,mask] - 
                 modelz*np.polyval(poly_coeff, alf_data[0,mask]) ) /alf_data[2,mask]

        if retval == 'resid':
            return resid

        chi2 = (resid**2).sum()

        if retval == 'lnp':        
            print(fit_params, chi2)
            if not np.isfinite(chi2):
                return -np.inf
            else:
                return -0.5*chi2

        return chi2


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
    
    