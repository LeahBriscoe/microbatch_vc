auc_fix = False
for n in names:
    print(n)
    for method in methods:
        print(method)
        if auc_fix:
            print(1-np.mean(all_methods_auc_stats[method][n]['auc_all']))
        else: 
            print(np.mean(all_methods_auc_stats[method][n]['auc_all']))



for name in names:    
    text_size = 10
    title = 'ROC curves'
    x_label = 'False positive rate'
    y_label = 'True positive rate'



    nplots=1
    plot_par_factor=2.26
    runs_cv_folds = 10
    plot_alpha=0.2
    fig, ax = plt.subplots(nplots, sharex=True, sharey=True)
    import seaborn as sns
    sns.set()
    current_palette = sns.color_palette()
    #sns.palplot(current_palette)
    palette = sns.color_palette("hls", len(methods))
    plot_color = palette.as_hex()
    plot_marker = ['None','None']
    plot_ls = ['-','--']
    plot_lw = 2
    plot_title = "test"

    for i in range(0,len(methods)):
        #'Naive Bayes''Random Forest'
        if auc_fix:
            fpr_all = 1- np.array(all_methods_auc_stats[methods[i]][name]['fpr_all'])
            tpr_i = 1- np.array(all_methods_auc_stats[methods[i]][name]['tpr_i'])
            tpr_all = 1-np.array(all_methods_auc_stats[methods[i]][name]['tpr_all'])
        else:
            fpr_all = all_methods_auc_stats[methods[i]][name]['fpr_all']
            tpr_i = all_methods_auc_stats[methods[i]][name]['tpr_i']
            tpr_all = all_methods_auc_stats[methods[i]][name]['tpr_all']

        ax.fill_between(fpr_all, tpr_all-np.std(tpr_i, axis=0)*plot_par_factor/np.sqrt(runs_cv_folds), tpr_all+np.std(tpr_i, axis=0)*plot_par_factor/np.sqrt(runs_cv_folds), color=plot_color[i], lw=0, alpha=plot_alpha)
        ax.plot(fpr_all, tpr_all, color=plot_color[i], ls='-', lw=plot_lw, marker=plot_marker[0])

    # fig.subplots_adjust(hspace=0)
    ax.set_xlabel(x_label, size=text_size)
    ax.tick_params(labelsize=text_size, axis='x')
    ax.set_ylabel(y_label, size=text_size)
    ax.tick_params(labelsize=text_size, axis='y')
    #ax.text(plot_title,plot_title,plot_title, va='center', ha='center', size=text_size+2)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    # ax.set_yticklabels(ax.get_yticks()[:-1])
    # ax.set_title(title, size=text_size+2)
    #leg = ax.legend()
    #leg_col = [plt.Rectangle((0, 0), 1, 1, fc=s, linewidth=0) for s in plot_color] + [plt.Line2D([0,1], [0,1], c='k', ls='-', lw=2)] + [plt.Line2D([0,1], [0,1], c='k', ls='--', lw=2)]
    #leg = ax.legend(leg_col, prop={'size':text_size}, loc='center left', bbox_to_anchor=(1.02,1), numpoints=1)
    #leg.get_frame().set_alpha(0)

    class LegendObject(object):
        def __init__(self, facecolor='white', edgecolor='white', dashed=False):
            self.facecolor = facecolor
            self.edgecolor = edgecolor
            self.dashed = dashed

        def legend_artist(self, legend, orig_handle, fontsize, handlebox):
            x0, y0 = handlebox.xdescent, handlebox.ydescent
            width, height = handlebox.width, handlebox.height
            patch = mpatches.Rectangle(
                # create a rectangle that is filled with color
                [x0, y0], width, height, facecolor=self.facecolor,
                # and whose edges are the faded color
                edgecolor=self.edgecolor, lw=3)
            handlebox.add_artist(patch)

            # if we're creating the legend for a dashed line,
            # manually add the dash in to our rectangle
            if self.dashed:
                patch1 = mpatches.Rectangle(
                    [x0 + 2*width/5, y0], width/5, height, facecolor=self.edgecolor,
                    transform=handlebox.get_transform())
                handlebox.add_artist(patch1)

            return patch

    from matplotlib.colors import colorConverter as cc
    import matplotlib.patches as mpatches
    bg = np.array([1, 1, 1])  # background of the legend is white
    colors = plot_color
    # with alpha = .5, the faded color is the average of the background and color
    colors_faded = [(np.array(cc.to_rgb(color)) + bg) / 2.0 for color in colors]

    ax.legend(list(range(6)), methods[0:6],
               handler_map={
                   0: LegendObject(colors[0], colors_faded[0]),
                   1: LegendObject(colors[1], colors_faded[1]),
                   2: LegendObject(colors[2], colors_faded[2]),
                   3: LegendObject(colors[3], colors_faded[3]),
                   4: LegendObject(colors[4], colors_faded[4]),
                   5: LegendObject(colors[5], colors_faded[5])
                })
    ax.set_facecolor('white')
    #ax.title('AUC confidence interval plot')
    fig.tight_layout()
    #fig.set_axis_on()
    ax.grid(color = "grey")
    fig.show()


    fig.savefig(str(plot_folder) + '/AUC_' + str(name) + str('.pdf'))