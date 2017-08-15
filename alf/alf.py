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
    def __init__(self, verbose=True):
                
        if not bool(driver.is_setup):
            if verbose:
                print('Alf: CALL SETUP()')
        
            driver.setup_alf()
        
        self.nl = driver.get_nspec()
        self.npar = driver.get_npar() 
        self.default_params = driver.get_default_parameters(self.npar)
        self.wave = driver.get_grid_lam(self.nl)
        self.params = self.default_params*1
        self.get_model()
        
    def get_model(self, in_place=True, **kwargs):
        self.set_param(**kwargs)
        sp = driver.get_spec(self.nl, self.params, self.npar)
        if in_place:
            self.spec = sp
        else:
            return sp
            
    def set_param(self, **kwargs):
        for k in kwargs:
            if k not in ALF_PARAMS:
                print('Paramter "{0}" not valid.'.format(k))
            else:
                self.params[pidx[k]] = kwargs[k] 
    
    def show_params(self):
        for k in ALF_PARAMS:
            print('{0:>15s} = {1:8.2f}'.format(k, self.params[pidx[k]]))
            