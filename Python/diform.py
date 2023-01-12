import os
import pickle
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.special import comb   
from scipy.spatial import ConvexHull
from matplotlib import pyplot as plt
import seaborn as sns
from Snoopy.Statistics import POT_GPD, rolling_declustering


class _DirectIformABC(  ) :
    """Abstract class for direct IFORM. Subclass are associated to a given univariate method. The univariate approach is implemented in the ._fit() method.
    
    Implementation of the approach described in  https://www.researchgate.net/publication/365459533_Model-free_environmental_contours_in_higher_dimensions
    
    More or less the equivalent of the following matlab code:  https://github.com/edmackay/Direct-IFORM
    """

    def __init__(self, df , npoints):
        """Abstract class, this will never be called directly

        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing n columns will result in coutour in n dimension
        npoints : int
            Discretisation for all dimensions

        """

        self._df = df.copy()
        self.npoints = npoints
        
        
        self.dim_names = self._df.columns.values
        self.ndim = len(self.dim_names)

        #--- Normalisation        
        self.mean_dict = { k : np.mean( df.loc[ : , k  ] ) for k in self.dim_names  }
        self.std_dict = { k : np.std( df.loc[ : , k  ] ) for k in self.dim_names  }
        for d in self.dim_names : 
            self._df[ d + "_a"  ] =  (df[ d  ] - self.mean_dict[ d] ) / self.std_dict[d]
        self.dim_names_a = [n+"_a" for n in self.dim_names]
        
        
        #--- Calculation of direction vectors
        self._direction_vector = direction_vector(self.ndim , self.npoints).rename( columns = { i:n for i, n in enumerate(self.dim_names_a) } )
        self._n_angles = len(self._direction_vector)
        self._res_df = pd.DataFrame( index = self._direction_vector.index, dtype = float )
        
        #-- Initialisation of array containing results of fit for each direction
        self._fit_objects = None 
        

    @property
    def results_df(self):
        return pd.concat( [self._direction_vector , self.fit_df, self._res_df  ], axis = 1 )
    
    @property
    def fit_df(self) : 
        return pd.DataFrame( index=self._direction_vector.index )
        
    def _fit(self , se):
        """Method to override in subclass. 

        Should return an object that has the following method :
          - rp_to_x
          - plot_rp_fit
          - plot_rp_data
        """
        raise(NotImplementedError())
        
        
    def to_pickle(self , filename):
        """Store to pickle format

        Parameters
        ----------
        filename : str
            File name
        """
        with open(filename , "wb") as f :
            pickle.dump( self, f )

    @classmethod
    def from_pickle(cls , filename):
        """Load from pickle
        
        Parameters
        ----------
        filename : str
            File name

        Returns
        -------
        DirectIform
            The DirectIform object
        """
        with open(filename , "rb") as f :
            o = pickle.load( f )
        return o
        
        

    def _scale_back(self , points_df ) :
        """Scale back the data

        Parameters
        ----------
        points_df : pd.DataFrame
            The dataframe to scale back

        Returns
        -------
        pd.DataFrame
            Dataframe in original dimension
        """
        return pd.DataFrame( index=points_df.index, data = { k : self.std_dict[k] * points_df[k+"_a"] + self.mean_dict[k]  for k in self.dim_names } )

        
    def fit_projections(self) :
        """Performs the univariate fits for all the direction vectors. 
        
        Results is stored in "_fit_objects"
        """
        
        if self._fit_objects is None : 
            _array = self._df.loc[ : , self.dim_names_a ].values
            _dir_vector = self._direction_vector.loc[:, self.dim_names_a].values
            self._fit_objects = np.empty( (self._n_angles) , dtype = object )
            for i, r in tqdm(list(enumerate(_dir_vector)), desc = "Fitting projections") :
                proj  = _array[:,:].dot( r )
                se =  pd.Series( index = self._df.index.values , data = proj)
                self._fit_objects[i] = self._fit( se )
                

    def rp_projections(self, rp_list ):
        """Calculate return values for on each projection, using existing fitted model
        
        Parameters
        ----------
        rp_list : np.ndarray
            List of return period
        """
        
        res = np.zeros( ( self._n_angles , len(rp_list) ) )
        for i in tqdm( np.arange(self._n_angles), desc = "Evaluate RP at each projection" ) :
            res[i,:] = self._fit_objects[i].rp_to_x( rp_list )

        for i , rp in enumerate(rp_list) :
            self._res_df[rp] = res[:,i]


    def _extract_contour( self , rp ) :
        """Extract contour at RP, using Voronoi cells, in working space (normalised)

        Parameters
        ----------
        rp : float
            Return period

        Returns
        -------
        np.ndarray
            Points
        """
        from scipy.spatial import Voronoi
        reflection = 2 * self._direction_vector.values * self._res_df.loc[:,rp].values[: , None]
        reflection = np.concatenate( [ np.array( [ np.zeros( (self.ndim) ),] ) , reflection])
        vor = Voronoi( reflection )
        return vor.vertices[ vor.regions[ vor.point_region[0] ] ]
    
    
    def extract_contour( self , rp ):
        """Extract contour (scaled back).
        
        Parameters
        ----------
        rp : Float
            Return period

        Returns
        -------
        pd.DataFrame
            Points of the contour
        """
        if rp not in self._res_df.columns :
            self.rp_projections( np.array([rp]) )
        return self._scale_back( pd.DataFrame(  data = self._extract_contour(rp=rp) , columns = self.dim_names_a ) )
       
        
       
    @staticmethod
    def _variable_change(df, output_variables, transform_dict):
        """Add columns according to operation specified in transform_dict

        Parameters
        ----------
        df : pd.DataFrame
            The input/output dataframe (modified!)
        output_variables : list
            The list of desired columns to get.
        transform_dict : dict
            Dictionary of functions

        Returns
        -------
        df : pd.DataFrame
            The dataframe containing all the output variables
        """
        for v in output_variables : 
            if v not in df.columns :
                df[v] = transform_dict[v] (df)
        return df
        

    def projected_contour( self , variables, rp, final_variables = None, return_triangulation = False, transform_dict = {} ):
        """Return projection of the contour

        Parameters
        ----------
        variables : list(str)
            Space in which the convex contour is calculated
        rp : float
            Return period
        final_variables : TYPE, optional
            Space in which the contour is output. The default is set to variables.
        return_triangulation : bool, optional
            If True, the triangulation is returned, necessary to display 3D surface. The default is False.
        transform_dict : dict, optional
            Dictionary of function to ease variable change. The default is {}.

        Returns
        -------
        pd.DataFrame
            The contour projection
        """        
        
        if final_variables is None : 
            final_variables = variables
        
        nd_contour = self.extract_contour(rp)
        
        self._variable_change( nd_contour, variables, transform_dict  )
       
        
        conv_hull = ConvexHull( nd_contour.loc[ : , variables ].values)
        
        # Construct dataframe containing only the contour points. 
        contour_points = pd.DataFrame( columns = variables, data = conv_hull.points[conv_hull.vertices] )
        
        
        self._variable_change( contour_points, final_variables, transform_dict  )
        
        output_points = contour_points.loc[: , final_variables]

        if return_triangulation :
            # ConvexHull.simplices return position in the original point vector. The two following lines convert to position of contour points
            _conv = pd.Series( index = conv_hull.vertices, data = np.arange( len( conv_hull.vertices ) ) )
            simplices = _conv.loc[ conv_hull.simplices.flatten() ].values.reshape( conv_hull.simplices.shape )
            return output_points, simplices
        else : 
            return output_points



    def sliced_contour(self, slice_dims , slice_values , rp, final_variables = None, transform_dict = {} ) :
        """Cut the contour at a given level.

        Parameters
        ----------
        slice_dim : str
            The cut dimension.
        slice_value : float
            The value at which to cut
        rp : float
            Return period
        return_triangulation : bool, optional
            If True, the triangulation is returned. The default is False.

        Returns
        -------
        pd.DataFrame
            The intersection points. (and optionally the triangulation)
        """
        
            
        nd_contour = self.extract_contour( rp )
        
        _slice = nd_contour
        
        for slice_dim, slice_value in zip( slice_dims, slice_values ) : 
            _slice = slice_convhull_df( _slice , slice_dim = slice_dim, slice_value=slice_value )       
            
            
        if final_variables is not None : 
            return self._variable_change( _slice, final_variables, transform_dict  ).loc[:,final_variables]
        else : 
            return _slice
            
            
            

    #------------ Plotting functions
    def plot_slice_2d(self,  slice_dims , slice_values , rp ,final_variables = None, ax = None, transform_dict = {}, color = "b", **kwargs):
        
        if len(slice_dims) != self.ndim-2 : 
            raise(Exception("Not a 2D slice"))
            
        slice_cont = self.sliced_contour( slice_dims , slice_values , rp  )
        

        label = kwargs.pop( "label" , f"RP = {rp:} " +  " ".join( f"{d:} = {v:}" for d, v in zip(slice_dims, slice_values)) )
        
        ax = self.plot_2d_convex_hull(slice_cont, final_variables=final_variables, transform_dict=transform_dict, ax=ax, color = color, label=label,**kwargs)

        if final_variables is None : 
            ax.set(xlabel = slice_cont.columns[0] , ylabel = slice_cont.columns[1] )
        else: 
            ax.set(xlabel = final_variables[0] , ylabel = final_variables[1] )
        return ax

    
    
    def plot_angle_parameters(self , plane, values = None, ax=None, **kwargs) :

        if ax is None :
            fig, ax = plt.subplots()

        if values is None : 
            values = self._res_df.columns
        id_ = np.isclose( (self.results_df.loc[ : , plane ]**2).sum(axis = 1) , 1 )
        res_df_sub = self.results_df.loc[id_].copy()

        res_df_sub["Angle"] =  np.rad2deg(np.mod( np.arctan2( res_df_sub[ plane[1] ] , res_df_sub[ plane[0] ] ) , 2*np.pi))
        res_df_sub.set_index("Angle").sort_index().loc[:, values ].plot(ax=ax, **kwargs)
        ax.set(ylabel = "Projected value", xlabel = f"Angle ( atan2( {plane[1]:} , {plane[0]:}) )" )
        return ax

    
    def plot_projection_2d( self , variables , rp ,final_variables = None, ax = None, transform_dict = {} , backend = "matplotlib", **kwargs ) : 
        """Plot the 2D projection of the contour
        
        Parameters
        ----------
        variables : list(str)
            Space in which the convex contour is calculated
        rp : float
            Return period
        final_variables : TYPE, optional
            Space in which the contour is plotted. The default is set to variables.
        transform_dict : dict, optional
            Dictionary of function to ease variable change. The default is {}.
        backend : str, optional
            Library used to plot. The default is "matplotlib".
        **kwargs : any
            Argument passed to the plotting function.
        Returns
        -------
        ax : TYPE
            The figure
        """
        
        cont = self.projected_contour(variables, rp, final_variables=final_variables, transform_dict = transform_dict)
        x, y = final_variables
        
        if backend == "matplotlib" :
            if ax is None :
                fig, ax = plt.subplots()
            ax.plot( cont.loc[:,x] , cont.loc[:,y] , **kwargs )
            ax.set(xlabel = cont.columns[0] , ylabel = cont.columns[1] )
            return ax
        elif backend == "plotly" :
            from plotly import express as px
            ax = px.line( data_frame = cont , x = x, y = y )
            return ax
        
    
    
    def plot_projection_3d(self, variables,  rp , ax = None, final_variables = None, transform_dict = {}, backend = "matplotlib", **kwargs):
        """Plot the 3D projection of the contour
        
        Parameters
        ----------
        variables : list(str)
            Space in which the convex contour is calculated
        rp : float
            Return period
        final_variables : TYPE, optional
            Space in which the contour is plotted. The default is set to variables.
        transform_dict : dict, optional
            Dictionary of function to ease variable change. The default is {}.
        backend : str, optional
            Library used to plot. The default is "matplotlib".
        **kwargs : any
            Argument passed to the plotting function.
        Returns
        -------
        ax : TYPE
            The figure
        """
        
           
        cont, tri = self.projected_contour(variables, rp, final_variables=final_variables, transform_dict = transform_dict, return_triangulation = True)

        tri = orient_normals( cont.values, tri )

        if backend == "plotly" :
            import plotly.figure_factory as ff
            
            fig = ff.create_trisurf(x=cont.iloc[:,0], y=cont.iloc[:,1], z=cont.iloc[:,2],
                                 simplices = tri ,
                                 plot_edges = False, **kwargs)
            fig.layout.scene.xaxis.title='Hs'
            fig.layout.scene.yaxis.title='Tm'
            fig.layout.scene.zaxis.title='U10'
            fig.show()
        else :
            if ax is None: 
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')           
            ax.plot_trisurf(  cont.iloc[:,0], cont.iloc[:,1], cont.iloc[:,2] , triangles = tri, **kwargs )
            ax.set(xlabel = "Hs" , ylabel ="U10" , zlabel = "Tm")
        return ax

        
    def plot_data_2d( self, variables , transform_dict = {}, ax = None, backend = "matplotlib", **kwargs ):
        """Plot the 2D projection of the contour
        
        Parameters
        ----------
        variables : list(str)
            Space in which the data are plotted.
        transform_dict : dict, optional
            Dictionary of function to ease variable change. The default is {}.
        backend : str, optional
            Library used to plot. The default is "matplotlib".
        **kwargs : any
            Argument passed to the plotting function.
        Returns
        -------
        ax : TYPE
            The figure
        """
                
        df = self._variable_change(self._df.copy(), variables, transform_dict)
        x, y = variables
        
        if len(variables) > 2 :
            color = variables[2]
        else:
            color = None

        if backend == "plotly" : 
            import plotly.express as px
            fig = px.scatter( data_frame = df , x = x , y = y, color = color, **kwargs )
            if ax is not None :
                import plotly.graph_objects as go
                fig = go.Figure(data = fig.data + ax.data , layout = ax.layout)
            return fig
        else : 
            if ax is None :
                fig, ax = plt.subplots()
            sns.scatterplot(data = df, x = x  , y = y  , hue = color , ax = ax, **kwargs)
            return ax


    def plot_univariate(self, i_direction, ax=None, fit_kwargs={}, data_kwargs={}):
        if ax is None :
            fig, ax = plt.subplots()
            
        self._fit_objects[i_direction].plot_rp_fit(ax=ax, **fit_kwargs)
        self._fit_objects[i_direction].plot_rp_data(ax=ax, **data_kwargs)
        
        t = [ f"{c:} = {self._direction_vector.loc[ 12, c ]:.2f}" for c in self._direction_vector.columns ]
        ax.set(title = ";".join(t) )
        return ax
        
        
        
        


    @staticmethod
    def plot_2d_convex_hull(points_df, final_variables = None, ax = None, transform_dict = {}, **kwargs):
        if ax is None : 
            fig, ax = plt.subplots()
        chull = ConvexHull(points_df.values)
        
        if final_variables is not None: 
            final_points = _DirectIformABC._variable_change( points_df,  final_variables, transform_dict = transform_dict ).loc[:,final_variables ].values
        else : 
            final_points =  chull.points
            
        label = kwargs.pop("label" , None)
        for i, d in enumerate(chull.simplices) : 
            ax.plot( final_points[d,0], final_points[d,1], label = label if i == 0 else None, **kwargs )
        return ax



def orient_normals( vertices, faces ) :
    """Orient the normals of a convex hull towards the exterior
    """
    faces = np.array(faces)
    tris = vertices[faces]
    
    # Calculate the normals
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    centers = np.mean( tris - np.mean(vertices) , axis = 1 )
    id_ = np.diag( np.matmul( centers, np.transpose(n) )   ) > 0
    faces_oriented = np.array(faces)
    faces_oriented[ id_ , 2 ]  = faces[ id_ , 1 ]
    faces_oriented[ id_ , 1 ]  = faces[ id_ , 2 ]
    return faces_oriented


        
class DirectIform( _DirectIformABC ) :
    """Specialisation of DirectIform when univariate fits are performed using POT and GP fit
    """
        
    def __init__(self, df, npoints, duration, window_int = None, window = None, threshold_q = 0.9, pot_kwargs = {}):
        """Compute direct IFORM contour, using POT with GP as univariate fit.

        Parameters
        ----------
        df : pd.DataFrame
            The input data
        npoints : int
            Discretisation for all dimensions.
        duration : float
            Duration corresponding to the input data
        window_int : int, optional
            Minimum window used to decluster the data. The default is None.
        window : TYPE, optional
            Minimum window used to decluster the data, in df.index scale. The default is None.
        threshold_q : float, optional
            Quantile to use for the threshold. The default is 0.9.
            
            
        Example
        -------
        >>> diform_2d = DirectIform( df.loc[: , ["Hs**0.5" , "Tm" ] ], npoints = 11 , window_int= 48 , duration = len(df) / (24*365), threshold_q = 0.9  )
        >>> diform_2d.fit_projections( )
        >>> contour_df = diform_2d.extract_contour( rp = 25.0 )
        """
        
        _DirectIformABC.__init__(self , df = df , npoints = npoints )
        self.window_int = window_int
        self.window = window
        self.threshold_q = threshold_q
        self.duration = duration
        self._pot_kwargs = pot_kwargs
        

    def _fit(self , se ) :
        
        declust = rolling_declustering(se , window_int = self.window_int , window = self.window)

        threshold = np.quantile(declust , self.threshold_q )

        pot = POT_GPD( declust , duration = self.duration, threshold = threshold, **self._pot_kwargs  )

        pot._fit()

        return pot


    @property
    def fit_df(self) : 
        return pd.DataFrame( index=self._direction_vector.index , data = { "THRESHOLD" : [ pot.threshold for pot in self._fit_objects], 
                                                                           "SCALE" : [ pot.scale for pot in self._fit_objects], 
                                                                           "SHAPE" : [ pot.shape for pot in self._fit_objects], 
                                                                           "NNLF" : [ pot.nnlf for pot in self._fit_objects], 
                                                                           "KS-PVALUE" : [ pot.ks for pot in self._fit_objects], 
                                                                         } )
    


def direction_vector( ndim, npoints, mirror = "add_negative" ) : 
    """Compute directtion vector to use for direct IFORM (or direct sampling) calculations
    """
    dims = list(range(ndim))
    df = pd.DataFrame( index = pd.MultiIndex.from_product( [ np.linspace(0,1, npoints) for n in range(ndim) ], names = dims ) , columns = ["sum"])
    df["sum"] = df.reset_index().loc[:, dims].sum(axis = 1).values
    res = df.loc[ np.isclose( df["sum"] , 1.0 ) , : ].reset_index().loc[ :,dims ]

    n_expected = comb(npoints+ndim-2,ndim-1, exact = True)
    
    if len( res ) != n_expected : 
        raise(Exception( "Problem in finding the correct number of point on the sphere of radius = 1" ))
      
    # Add negative side if required
    if type(mirror) == str : 
        mirror = [mirror for i in range(ndim)]
    
    for idim in range(ndim) : 
        if mirror[idim] == "add_negative" : 
            dup = res.loc[ res.loc[ : , idim ] > 0 ].copy()
            dup.loc[: , idim ] *=-1
            res = pd.concat( [res, dup ]  )
    res.reset_index(inplace = True, drop = True)
            
    # Normalize
    res.values[:,:] /= ((res**2).sum(axis = 1)**0.5).values[:,None]
    
    return res

        
def slice_convhull( conv_hull, cut_dim, slice_value ):
    """Calculate intersection point of a convex hull

    Parameters
    ----------
    conv_hull : ConvexHull
        The hyper-surface to cut
    cut_dim : int
        The dimension in which the cut is performed
    slice_value : float
        Value at which to slice
        
    Returns
    -------
    intersection_point : np.ndarray
        The intersection point
    """
        
    # Find cutted simplices
    x_min = conv_hull.points[ conv_hull.simplices , cut_dim  ].min(axis = 1)
    x_max = conv_hull.points[ conv_hull.simplices , cut_dim  ].max(axis = 1)
    cutted_simplices = conv_hull.simplices[ np.where( (x_min < slice_value) &  (x_max >= slice_value) ) ]
    
    # Find edges of cutted simplices
    permutation = list(itertools.combinations( np.arange(conv_hull.ndim), 2))
    edges = np.concatenate(  [ cutted_simplices[: , permutation[i]] for i in range(len( permutation )) ] )
    
    # Find cutted edges of cutted simplices
    x_min_e = conv_hull.points[ edges, cut_dim ].min(axis = 1)
    x_max_e = conv_hull.points[ edges, cut_dim ].max(axis = 1)
    cutted_edges = edges[ np.where( (x_min_e < slice_value) &  (x_max_e >= slice_value) ) ]


    #Interpolate the intersection points    
    alpha =  (conv_hull.points[ cutted_edges[:,1], cut_dim ] - slice_value) / np.diff( conv_hull.points[ cutted_edges, cut_dim ], axis = 1)[:,0]
    intersection_point = np.zeros( (len(cutted_edges) , conv_hull.ndim-1) )
    for idim in range(conv_hull.ndim-1):
        intersection_point[:, idim] = alpha * conv_hull.points[ cutted_edges[:,0], idim ]   + (1-alpha)*conv_hull.points[ cutted_edges[:,1], idim ] 
        
    intersection_point = np.unique( intersection_point, axis = 0 )
    
    return intersection_point


def slice_convhull_df( points_df, slice_dim, slice_value ):
    """Same as slice_convhull, but with dataframe as input and output
    """
    conv_hull = ConvexHull( points_df )
    points = slice_convhull( conv_hull , cut_dim = points_df.columns.get_loc(slice_dim), slice_value=slice_value )       
    contour_slice = pd.DataFrame( columns = points_df.columns.drop(slice_dim).values , data = points )
    return contour_slice




    

    
