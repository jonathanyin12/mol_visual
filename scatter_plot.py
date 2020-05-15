from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.io import output_notebook
from bokeh.palettes import Category20
import numpy as np

def plot_mols(x, y, labels, imgs, groups, s=8, alpha = 0.8, x_label='', y_label='', title='', size=(800,600)):  
    TOOLTIPS = """
        <div>
            <div>
                <span style="font-size: 12px; font-weight: bold;">@name</span>
            </div>
            <div>
                <img src="@imgs" height="100" alt="@imgs" border="1.5"></img>
            </div>
        </div>
    """
    p = figure(plot_width=size[0], 
               plot_height=size[1],  
               tooltips=TOOLTIPS, 
               title=title,
               x_axis_label=x_label, 
               y_axis_label=y_label)

    colors = Category20[20][:len(np.unique(groups))]
    
    points = ColumnDataSource(data=dict(
       x=x,
       y=y,
       name=labels,
       imgs=imgs,
       color=[colors[x] for x in groups]
    ))
        
    p.circle('x', 'y',  fill_color='color',line_color=None, size=s, alpha = alpha, source=points)

    output_notebook()
    output_file("toolbar.html")
    return p
