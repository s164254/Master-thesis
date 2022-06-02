import matplotlib.pyplot as plt

def dataframe_plot(df, plot_func, title, axis_setup_func=None, title_font_size=4, tick_font_size=4, block=False, fig_filename=None):
    ax = plot_func(df)
    
    if axis_setup_func:
        axis_setup_func(ax)

    plt.xticks(fontsize=tick_font_size)
    plt.yticks(fontsize=tick_font_size)
    plt.title(label=title, fontsize=title_font_size)
    do_block = block
    if fig_filename:
        do_block = False
    plt.show(block=do_block)
    if fig_filename:
        plt.savefig(fig_filename)
