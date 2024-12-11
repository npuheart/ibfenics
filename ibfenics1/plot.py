import matplotlib.pyplot as plt

# plot_multiple_lines(time, data, labels,types, title='porous pressure', xlabel='time step', ylabel='p')
def plot_multiple_lines_1(
    time,
    data,
    labels=None,
    types=None,
    colors=None,
    linestyles=None,
    title="",
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    xticks=None,
    yticks=None,
    ncol=1,
    legends=None,
    ylog=False,
):
    """
    Plot multiple lines from a list of lists.
    
    Parameters:
        data (list of lists): List of data points for each line.
        labels (list): List of labels for each line. Default is None.
        title (str): Title of the plot. Default is None.
        xlabel (str): Label for the x-axis. Default is None.
        ylabel (str): Label for the y-axis. Default is None.
    """
    basic_colors = [
        "#FC8002",
        "#369F2D",
        "#B9181A",
        "#1663A9",
        "#614099",
        "#369F2D",
        "#B9181A",
        "#1663A9",
    ]
    # basic_linestyles = [":",":",":","-","-","-","--","--","--"]
    basic_linestyles = [":", "-.", "--", "-", "-", "-", "--", "--", "--"]
    if not labels:
        labels = [f"Line {i+1}" for i in range(len(data))]

    if not types:
        types = ["lines" for i in range(len(data))]

    if not linestyles:
        linestyles = [
            basic_linestyles[i % len(basic_linestyles)] for i in range(len(data))
        ]

    if not colors:
        print("colors is None")
        print(basic_colors)
        print(len(data))
        colors = [basic_colors[i] for i in range(len(data))]
    
    # TODO: labels, types, 的列表长度要和data一样长
    for i in range(len(data)):
        if types[i] == "lines":
            plt.plot(
                time[i],
                data[i],
                # marker="o",
                color=colors[i],
                linestyle=linestyles[i],
            )
        elif types[i] == "dots":
            plt.scatter(time[i], data[i], color=colors[i], label=labels[i])
        else:
            print("画图的 types 参数错误")

    # 添加图例并指定图例在图外的位置, 并调整图形的边距以适应图例
    # legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.subplots_adjust(right=0.75)
    # legend = plt.legend()
    # legend = plt.legend(ncol=ncol)
    if legends:
        plt.legend(legends,loc='upper right')

    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    if xticks:
        plt.xticks(ticks=xticks)

    if yticks:
        plt.yticks(ticks=yticks)

    # plt.xticks(ticks=[i*0.5 for i in range(5)])
    # plt.xticks(ticks=[i*0.5 for i in range(5)])

    # 使用 bbox_extra_artists 包含图例
    # bbox_extra_artists=(legend,), bbox_inches='tight')

    ax = plt.gca()
    # 设置y轴为对数坐标
    if ylog:
        ax.set_yscale("log")

    # ax.set_aspect(1.1)
    # ax.set_aspect(1.0/(ylim[1]-ylim[0]))
    

    # 保存图像
    plt.savefig("figure.png",  bbox_inches='tight', dpi=300)
    plt.savefig("figure.svg",  bbox_inches='tight')
    plt.close()
