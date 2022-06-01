import matplotlib.pyplot as plt

def dataframe_plot(df, plot_func, title, axis_setup_func=None, title_font_size=4, tick_font_size=4, block=False):
    ax = plot_func(df)
    
    if axis_setup_func:
        axis_setup_func(ax)

    plt.xticks(fontsize=tick_font_size)
    plt.yticks(fontsize=tick_font_size)
    plt.title(label=title, fontsize=title_font_size)
    plt.show(block=block)
