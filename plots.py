import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

def frequency_plot(x_label, x_data):
    plt.figure(figsize=(5.5, 5.5))

    sns.countplot(x=x_label, data=x_data, edgecolor='black')

    plt.xlabel(f'Bioactivity {x_label}', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency', fontsize=14, fontweight='bold')

    plt.savefig(f'plots/frequency_plot_bioactivity_{x_label}.png')

def scatter_plot(x_label, y_label, data, hue, size):
    plt.figure(figsize=(5.5, 5.5))

    sns.scatterplot(x=x_label, y=y_label, data=data, hue=hue, size=size, edgecolor='black', alpha=0.7)

    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.savefig(f'plots/scatter_plot_{x_label}_vs_{y_label}.png')  

def box_plot(x_label,y_label,data):
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x = x_label, y = y_label, data = data)

    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')

    plt.savefig(f'plots/box_plot_{y_label}.png')

def regplot_plot(x_label, y_label, Y_test, Y_pred):
    plt.figure(figsize=(5.5, 5.5))
    sns.set(color_codes=True)
    sns.set_style("white")

    ax = sns.regplot(x=Y_test, y=Y_pred, scatter_kws={'alpha': 0.4})
    ax.set_xlabel(x_label, fontsize='large', fontweight='bold')
    ax.set_ylabel(y_label, fontsize='large', fontweight='bold')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    plt.savefig(f'plots/regplot_plot_{y_label.replace(" ", "")}.png')