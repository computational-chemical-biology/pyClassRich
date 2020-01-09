import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.style.use('ggplot')

def plot_volcano(vdf,  palette, mode='show', filename=None):
    ax = sns.scatterplot(x="fold_change", y="p-value",
                         hue="Class", data=vdf,
                         palette=palette)
    ax.set(xlabel='log 2 fold change', ylabel='-log10 p-value')
    if mode=='show':
        plt.show()
    elif mode=='pdf':
        with PdfPages(filename) as pdf:
            pdf.savefig()
            plt.close()


