import numpy as np
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import pyplot as plt
from scipy.optimize import minimize


def rolling_declustering( se , window=None, window_int=None ) : 
    """Return declustered events
   
    Parameters
    ----------
    se : pd.Series
        Time series (time is index)
    window : float, optional
        window used to decluster the data. The default is None.
    window_int : int, optional
        window used to decluster the data, in number of time step. The default is None.

    Returns
    -------
    pd.Series
        The declustered sample
    """
  
    if window_int is not None : 
        _se = se.reset_index(drop=True)
        se_tmp = _se.rolling( window = window_int, min_periods=1, center=True, axis=0, closed = 'neither' ).max()   
        se_tmp = se_tmp.loc[ se_tmp == _se ]
        se_tmp.loc[ np.concatenate( [[True], np.diff(se_tmp.index) > window_int / 2 ]) ]
        se_tmp.index = se.index[se_tmp.index.values]
    else :
        se_tmp = se.rolling( window = window, min_periods=1, center=True, axis=0, closed = 'neither' ).max()   
        se_tmp = se_tmp.loc[ se_tmp == se ]
        se_tmp = se_tmp.loc[ np.concatenate( [[True], np.diff(se_tmp.index) > window / 2] ) ]
    return se_tmp





def probN( n, variant = (0.0, 0.0), alphap = None, betap = None):
    """Return exceedance probability of the ranked data.

    (inverse of scipy.stats.mstats.mquantiles)


    Parameters
    ----------
    n : int
        Size of the data vector

    variant : float, tuple or str. 
        Variant for plotting positions parameter. The default is i / (n+1). 

    Returns
    -------
    np.ndarray
        Exceedance probability of the ranked data.

    Note
    ----
    If variant is a tuple (alphap , betap):
    
        Typical values of (alphap,betap) are:
            - (0,1)    : ``p(k) = k/n`` : linear interpolation of cdf
              (**R** type 4)
            - (.5,.5)  : ``p(k) = (k - 1/2.)/n`` : piecewise linear function
              (**R** type 5)
            - (0,0)    : ``p(k) = k/(n+1)`` :
              (**R** type 6)
            - (1,1)    : ``p(k) = (k-1)/(n-1)``: p(k) = mode[F(x[k])].
              (**R** type 7, **R** default)
            - (1/3,1/3): ``p(k) = (k-1/3)/(n+1/3)``: Then p(k) ~ median[F(x[k])].
              The resulting quantile estimates are approximately median-unbiased
              regardless of the distribution of x.
              (**R** type 8)
            - (3/8,3/8): ``p(k) = (k-3/8)/(n+1/4)``: Blom.
              The resulting quantile estimates are approximately unbiased
              if x is normally distributed
              (**R** type 9)
            - (.4,.4)  : approximately quantile unbiased (Cunnane)
            - (.35,.35): APL, used with PWM
            
    if variant is a float: 
            p = (i+(a-1)/2) / (N+a)        
    """
       
    k = np.arange(1, n+1 , 1)    
    
    if isinstance( variant , tuple ) or isinstance( variant , list ) : 
        alphap , betap = variant        
        return 1 - (k - alphap)/(n + 1 - alphap - betap)   
    elif isinstance( variant, float) :
        return 1 - ( (k + 0.5*(variant-1) ) / (n + variant)  )
    elif variant == "median" : 
        b = stats.beta( k , n - k + 1 )
        return 1 - b.median()
    else : 
        raise( "Unknown variant : {variant:}" )


def probN_ci( n, alpha = 0.05):
    """Exceedance probability CI of a sorted sample of size n. 
    
    Parameters
    ----------
    n : int
        Number of samples
    alpha : float, optional
        1 - Size of the confidence interval. The default is 0.05.
        
    Returns
    -------
    ci_l : np.ndarray
        Lower bound of the CI
    ci_u : np.ndarray
        Upper bound of the CI
    """
    m = np.arange( 1, n+1 )
    ci_l = stats.beta( m , 1 + n - m ).ppf( alpha/2 )
    ci_u = stats.beta( m , 1 + n - m ).ppf(  1-alpha/2  )
    return 1-ci_u, 1-ci_l


class POT():

    def __init__(self, sample, duration, threshold, variant = (0.,0.) ):
        """Peak over Threshold method to calcualte return values. 
        
        Uses only empirical quantiles, for generalized pareto fit, see POT_GPD class. 

        Parameters
        ----------
        sample : np.ndarray
            Sample of independant observation
        duration : float
            Duration corresponding to the sample
        threshold : float
            Threshold
        """
        
        self.sample = sample
        
        self.duration = duration
        
        self.threshold = threshold
        
        self.extremes = np.sort( sample[sample >= threshold] )
        self.exceedances = self.extremes - self.threshold
        
        self.f = len(self.extremes) / self.duration
        
        #Which variant for the empirical quantile calculation
        self._variant = variant
        
        #Interpolator, not always needed ==> lazy evaluated
        self._x_to_rp_empirical = None
        self._rp_to_x_empirical = None
        
        #POT object using resampled data
        self._bootstrap_objects = []
        
        
    def x_to_rp_empirical(self, x):
        """Return period from return value, using interpolation between data.

        Parameters
        ----------
        x : float
            Return value

        Returns
        -------
        float
            Return period
        """
        if self._x_to_rp_empirical is None : 
            self._build_interp_x_to_rp()
        return self._x_to_rp_empirical(x)
    
    def rp_to_x_empirical(self, x):
        """Return value from return period

        Parameters
        ----------
        rp : float
            return period

        Returns
        -------
        float
            return value
        """
        if self._rp_to_x_empirical is None : 
            self._build_interp_rp_to_x()
        return self._rp_to_x_empirical(x)    
            

    def _build_interp_rp_to_x(self):
        # Build interpolator
        self._rp_to_x_empirical = InterpolatedUnivariateSpline( self.empirical_rp() , self.extremes, ext = "raise", k = 1 )

    def _build_interp_x_to_rp(self):
        # Build interpolator
        self._x_to_rp = InterpolatedUnivariateSpline( self.extremes , self.empirical_rp(), ext = "raise", k = 1)


    @classmethod
    def FromTimeSeries( cls, se, duration, threshold = None, threshold_q = None, window = None , window_int = None, **kwargs ):
        """Create POT analysis using time series as input
        
        Parameters
        ----------
        se : pd.Series
            The time signal.
        duration : float
            Duration associated to the time-series.
        threshold : float
            Threshold.
        treshold_q : float
                Threshold.
        window : float, optional
            window used to decluster the data. The default is None.
        window_int : int, optional
            window used to decluster the data, in number of time step. The default is None.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        The POT analysis class
        """
        sample = rolling_declustering( se, window = window, window_int = window_int )
        
        if threshold_q is not None :
            threshold = np.quantile( sample, threshold_q, method = "weibull") # Weibull corresponds to "i / (n+1)" method
        return cls( sample.values, duration, threshold = threshold, **kwargs )
    
        
    def empirical_rp( self  ):
        """Return empirlcal return period of events above threshold (sorted).

        Parameters
        ----------
        variant : (float,float) or str, optional
            DESCRIPTION. The default is (0.0, 0.0), which corresponds to i/(n+1)

        Returns
        -------
        np.ndarray
            Return period of events above threshold (sorted).
        """
        return 1 / ( self.f * probN( len(self.extremes), variant = self._variant ) )


    def plot_rp_data(self, ax = None, variant = (0.0 , 0.0), marker = "+", linestyle = "", **kwargs):
        """
        

        Parameters
        ----------
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        variant : TYPE, optional
            DESCRIPTION. The default is (0.0 , 0.0).
        marker : TYPE, optional
            DESCRIPTION. The default is "+".
        linestyle : TYPE, optional
            DESCRIPTION. The default is "".

        Returns
        -------
        ax : TYPE
            DESCRIPTION.

        """
        
        if ax is None :
            fig, ax = plt.subplots()

        ax.plot( self.empirical_rp(),  self.extremes, marker = marker , linestyle=linestyle, **kwargs)
        ax.set_xscale("log")
        return ax
    
        

# Faster than scipy.stats implementation. (no 'loc' parameter in our case), TODO : pass to c++ for efficiency ?
_xmax_log = np.log(np.finfo(float).max)
def nnlf_genpareto( shape_scale, data ):
    shape, scale = shape_scale
    good = np.where( data * shape  / scale > -1 + 1e-12  )[0]
    ngood = len(good)
    n = len(data)
    nbad = n - ngood
    if nbad > 0 : 
        penalty = _xmax_log * nbad
        return penalty + ngood*np.log(scale) + ( 1. + 1. / shape)  * np.sum( np.log1p( shape * data[good] / scale) )
    else :
        return n*np.log(scale) + ( 1. + 1. / shape)  * np.sum( np.log1p( shape * data / scale) )
    



class POT_GPD( POT ):
    
    def __init__(self, sample, duration , threshold , shape_bounds = (-np.inf, 1.0), scale_bounds = (1e-12 , np.inf), solver = "minimize_mle", fit_kwargs = {} ):
        """Peak over threshold, extremes are fitted with Generalize Pareto Distribution
        
        Parameters
        ----------
        sample : np.ndarray
            Sample of independant observation
        duration : float
            Duration corresponding to the sample
        threshold : float
            Threshold
        shape_bounds : tuple, optional
            Bounds for the shape parameter
        solver : str, optional
            Which library to use for minimizing the likelyhood
        fit_kwargs : dict, optional
            Argument pass to scipy.stats.rv_continous.fit or scipy.optimize.minimize
        """
        POT.__init__(self , sample, duration , threshold)
        
        self._solver = solver

        # Update default minimizer settings        
        self._fit_kwargs = {}
        if self._solver == "minimize_mle":
            self._fit_kwargs = {"method" : "nelder-mead"}
        self._fit_kwargs.update( fit_kwargs )
        
        self.shape_bounds = shape_bounds
        self.scale_bounds = scale_bounds

        # Lazy evaluation
        self._gpd = None
        

    
    def _fit(self):
        # Start coefficients, from method of moments        
        mu = np.mean(self.exceedances)
        s = np.std(self.exceedances)
        _shape = 0.5 * (1-( mu/s )**2  )
        _scale = mu * ( 1-_shape)
        
        if self._solver == "minimize_mle" :  # Generallyt faster than stats.genpareto.fit
            
            res = minimize( nnlf_genpareto, x0 = [np.clip( _shape, *self.shape_bounds) , np.clip(_scale, *self.scale_bounds) ], args = (self.exceedances,), bounds = ( self.shape_bounds, self.scale_bounds ),  **self._fit_kwargs )
            self._shape , self._scale = res.x
        elif self._solver == "genpareto.fit" : 
            self._shape , _ , self._scale = stats.genpareto.fit( self.exceedances, _shape, scale = _scale, floc=0, **self._fit_kwargs)
        elif self._solver == "mom" : 
            self._shape , self._scale = _shape, _scale
        else : 
            raise(Exception())

        self._gpd = stats.genpareto(self._shape , 0.0 , self._scale)
        
        
        
    @property
    def gpd(self):
        if self._gpd is None:
            self._fit()
        return self._gpd

    @property
    def shape(self):
        return self.gpd.args[0]

    @property
    def scale(self):
        return self.gpd.args[2]
    
    
    @property
    def nnlf(self):
        return nnlf_genpareto( [self._shape , self.scale] , self.exceedances )
    
    @property
    def ks(self):
        return stats.kstest(self.exceedances , self.gpd.cdf ).pvalue
    
    def x_to_rp( self, x ) :
        """Calculate return period from return value

        Parameters
        ----------
        x : float or np.ndarray
            Return value

        Returns
        -------
        float or np.ndarray
            return period
        """
        return  1 /  (self.f * ( self.gpd.sf( x - self.threshold  ))  )
        
    def rp_to_x(self , rp):
        """Provide return value at RP
        
        Parameters
        ----------
        rp : float or array
            Return period.

        Returns
        -------
        float or np.ndarray
             Return value
        """
        return self.threshold + self.gpd.ppf( 1. - ( 1 / (rp * self.f ) ))
    
    

    def bootstrap(self, n = 1000):
        """Bootstrap the POT analysis

        Parameters
        ----------
        n : int, optional
            Number of re-sample. The default is 1000.
        """
        for i in range(n) : 
            self._bootstrap_objects.append( self.__class__( np.random.choice( self.sample , size = len(self.sample)  ) ,
                                                          duration = self.duration,
                                                          threshold = self.threshold,
                                                          fit_kwargs = self._fit_kwargs,
                                                          shape_bounds = self.shape_bounds,
                                                          scale_bounds = self.scale_bounds ))

   
    def rp_to_xci(self, rp, alpha_ci, ci_type = "bootstrap", n_bootstrap = 1000):
        """Return lower and upper bound of the confidence interval.

        Parameters
        ----------
        rp : float or array
            Return period.
        alpha_ci : float
            Centered confidence interval.
        ci_type : str, optional
            How the CI is evaluated. The default is "bootstrap".
        n_bootstrap : int, optional
            Number of re-sample for the bootstrap. The default is 1000.

        Returns
        -------
        x_low : float or array
            Lower bound of the confidence interval.
        x_up : float or array
            Upper bound of the confidence interval.
        """
        
        if len(self._bootstrap_objects) < n_bootstrap :
            self.bootstrap( n_bootstrap - len(self._bootstrap_objects) )

        if ci_type == "bootstrap":
            v = [ b.rp_to_x( np.array(rp)) for b in self._bootstrap_objects ]
            x_low = np.quantile(v , alpha_ci/2 , axis = 0, method = "weibull")
            x_up = np.quantile(v , 1-alpha_ci/2 , axis = 0, method = "weibull")
            return x_low, x_up
        else : 
            raise(Exception( f"ci_type {ci_type:} is not known" ))
    
    
    def plot_rp_fit(self, rp_range=None, ax=None, **kwargs):
        """Plot return value against return period.
        
        Parameters
        ----------
        rp_range : np.ndarray or None, optional
            Range of RP to plot. The default is None.
        ax : plt.Axis, optional
            The figure. The default is None.

        Returns
        -------
        plt.Axis
            The figure
        """
        
        if ax is None :
            fig, ax= plt.subplots()
            
        if rp_range is None : 
            _x = self.empirical_rp()
            rp_range = np.logspace(  np.log10( np.min( _x ) ) , np.log10(np.max( _x ))*1.5   , 200 )
        
        ax.plot( rp_range , self.rp_to_x( rp_range ), **kwargs )
        ax.set_xscale("log")
        ax.set( xlabel = "Return period" )
        return ax
    
    
    def plot_rp_ci( self, alpha_ci, rp_range = None, ax=None, ci_type = "bootstrap", alpha = 0.4, **kwargs) : 
        """Plot return value CI against return period.
        
        Parameters
        ----------
        rp_range : np.ndarray or None, optional
            Range of RP to plot. The default is None.
        alpha_ci : float
            Centered confidence interval.
        ci_type : str, optional
            How the CI is evaluated. The default is "bootstrap".
        ax : plt.Axis, optional
            The figure. The default is None.

        Returns
        -------
        plt.Axis
            The figure
        """
        if ax is None :
            fig, ax= plt.subplots()
            
        if rp_range is None : 
            _x = self.empirical_rp()
            rp_range = np.logspace(  np.log10( np.min( _x ) ) , np.log10(np.max( _x ))*1.5   , 200 )
            
        x_low, x_high = self.rp_to_xci( rp = rp_range, alpha_ci=alpha_ci , ci_type = ci_type )
        
        ax.fill_between( rp_range, x_low, x_high, alpha = alpha, **kwargs)
            
        return ax


if __name__ == "__main__" : 

    import pandas as pd
    from scipy.io import loadmat
    a = loadmat( "../MATLAB/Example/Data_Celtic_Sea.mat", struct_as_record = False )
    names = ["Hs" , "Tm" , "U10" , "dir_rel"]
    df = pd.DataFrame( index = a["time"][:,0] , data = { k:a[k][:,0] for k in names  } )
    
    data = df.Hs
    
    duration = len(data) / (365*24)
   
    pot = POT_GPD.FromTimeSeries( se = data, duration=duration, threshold_q = 0.9, window_int=48, solver = "minimize_mle" )
    pot.bootstrap(n=1000)  
    fig, ax = plt.subplots()
    pot.plot_rp_data(ax=ax, label = "Hs data", color = "blue")
    pot.plot_rp_fit(ax=ax, color = "blue")
    ax.legend()
    
    
    pot.plot_rp_ci(alpha_ci = 0.05, ax=ax, color = "lightblue")

    
    pot.rp_to_x_empirical( 10 )
    pot.rp_to_x( 10 )
    
    
