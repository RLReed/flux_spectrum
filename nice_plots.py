def init_nice_plots() :
    """ Some things to make plots that Roberts seems to like.
    """
    from matplotlib import rc
    rc('font',**{'family':'serif'})
    from matplotlib import rcParams
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    rcParams['xtick.labelsize'] = 18
    rcParams['ytick.labelsize'] = 18
    rcParams['lines.linewidth'] = 1
    rcParams['axes.labelsize'] = 20
    rcParams.update({'figure.autolayout': True}) 