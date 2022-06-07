import matplotlib.pyplot as plt
from pyparsing import Forward
import pylab as plot

#def dataframe_plot(df, plot_func, title, axis_setup_func=None, title_font_size=4, tick_font_size=4, block=False, fig_filename=None):
#    _dataframe_plot(df, plot_func, title, axis_setup_func, title_font_size, tick_font_size, block, fig_filename)

#def dataframe_plot_noxticks(df, plot_func, title, axis_setup_func=None, title_font_size=4, tick_font_size=4, block=False, fig_filename=None):
#    _dataframe_plot(df, plot_func, title, axis_setup_func, title_font_size, tick_font_size, block, fig_filename, xticks_on=False)


def dataframe_plot(df,
                   plot_func,
                   title,
                   axis_setup_func=None,
                   plot_setup_func=None,
                   title_font_size=8,
                   tick_font_size=8,
                   block=False,
                   fig_filename=None,
                   xlabel='',
                   ylabel='',
                   ylim=None,
                   xlim=None):
    params = {'legend.fontsize': 6} #, 'legend.handlelength': 2}
    plot.rcParams.update(params)

    ax = plot_func(df)

    if ylim:
        ax.set_ylim(*ylim)

    if xlim:
        ax.set_xlim(*xlim)

    # if legend:
    #     plt.legend(fontsize="medium") # using a named size
    #plt.rc('legend',fontsize='small') # using a named size

    plt.xticks(fontsize=tick_font_size)
    plt.yticks(fontsize=tick_font_size)
    plt.title(label=title, fontsize=title_font_size)

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    if axis_setup_func:
        axis_setup_func(ax)

    if plot_setup_func:
        plot_setup_func(plt)

    do_block = block
    if fig_filename:
        do_block = False
    plt.show(block=do_block)
    if fig_filename:
        print(fig_filename)
        plt.subplots_adjust(bottom=0.15)
        plt.rcParams['savefig.dpi'] = 300
        #plt.rcParams["figure.figsize"] = (10,6) # width, height
        plt.savefig(fig_filename)
