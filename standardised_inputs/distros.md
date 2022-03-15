## Distributions

```python
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.signal import savgol_filter
step = 0.1
bin_values = np.arange(start=-10, stop=10+step, step=step)
fig = go.Figure()
colors = ('#F8766D','#00BA38','#619CFF')
for i, filename in enumerate(('star.standarised.csv', 'O2567.standarised.csv', 'S1342.standarised.csv')):
    name = filename.split('.')[0]
    reference = pd.read_csv(filename, index_col=0)
    counts, bins = np.histogram(reference.experimental_ddG,
                                bins=bin_values)
    fig.add_trace(go.Scatter(x=0.5 * (bins[:-1] + bins[1:]), 
                             y=savgol_filter(counts/len(reference), 5, 3), # window size 11, polynomial order 3
                             name=name,
                             line=dict(width=3, color=colors[i]),
                             #opacity=0.1,
                             #fill='tozeroy' # fill down to xaxis
                            )
                 )
    fig.add_vline(x=reference.experimental_ddG.mean(),line_width=2, line_dash="dash", line_color=colors[i])
    fig.add_vline(x=reference.experimental_ddG.median(), line_width=2, line_dash="dot", line_color=colors[i])
    
fig.update_layout(title=f'Distribution of experimental ∆∆G in datasets<br>(Savitzky–Golay filtered with a window size of 10 bins, ie. {step*10} kcal/mol)',
                  paper_bgcolor='white',
                  plot_bgcolor='white', # or rgba(0,0,0,0)
                  xaxis={'title':'∆∆G [kcal/mol]', 'gridcolor': 'gainsboro', 'zerolinecolor': 'gray'},
                  yaxis={'title':f'Fraction of counts per {step} kcal/mol bin', 'gridcolor': 'gainsboro', 'zerolinecolor': 'gray'})
fig.write_image("distribution.png",scale=3)
fig.show()
```