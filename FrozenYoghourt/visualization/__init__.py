from FrozenYoghourt import *

def scatter(f, lb = 0, ub = 2*np.pi, data_points = 100, point_size = 2, title = '', 
            xlabel = '', ylabel = '', xscale = 'linear', yscale = 'linear', legend_list = []):
    
    """
    Draw scatter plot from the inputted function/s
    
    Parameters
    ----------
    f: function, list
        Individual function or list of functions to plot
    lb: float, optional
        Lower bound of the plot domain (default is 0)
    ub: float, optional
        Upper bound of the plot domain (default is 2pi)
    data_points: int, optional
        Number of scatter points (default is 100)
    point_size: int, optional
        Size of a scatter point (default is 2)
    title: str, optional
        Title of the plot (default is '')
    xlabel: str, optional
        Label of the x-axis (default is '')
    ylabel: str, optional
        Label of the y-axis (default is '')
    xscale: str, optional
        Scale of the x-axis (default is 'linear')
    yscale: str, optional
        Scale of the y-axis (default is 'linear')
    legend_list: list, optional
        List of names for the different plots (default is [])
        
    Returns
    -------
    fig: matplotlib.figure.Figure
        Figure of the plot
    ax: matplotlib.axes._subplots.AxesSubplot
        Axis of the plot
        
    """
    
    domain = np.linspace(lb, ub, data_points)
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_xscale(xscale); ax.set_yscale(yscale)

    if (type(f) == list) or (type(f) == tuple):
        codomain_list = [g(domain) for g in f]
        for index, codomain in enumerate(codomain_list):
            ax.scatter(domain, codomain, s = point_size)

        if legend_list == []:
            legend_list = [f'Function {i}' for i in range(1, len(f)+1)]
        ax.legend(legend_list)

    else:
        codomain = f(domain)
        ax.scatter(domain, codomain, s = point_size)
        ax.legend(legend_list)
        
    return fig, ax