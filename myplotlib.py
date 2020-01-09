import plotly
import plotly.graph_objs as go
from datetime import datetime
class PlotLib(object):
    """ get box_plot for a dict like data"""

    def __init__(self,data,):
        
        pass

    def get_str_datetime(self):
        return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def scatter(self, data, mode):
        pass

    def box(self,data, title,filename=None):
        if not filename:
            filename = self.get_str_datetime() + '_box.html'
        traces = []
        for key,value in data.items():
            trace = go.Box(
                y = value,
                name = key
            )
            traces.append(trace)
        layout = go.Layout(
            title = title,
            yaxis=dict(title='AUC')
        )
        fig = go.Figure(data=traces,layout=layout)
        plotly.offline.plot(fig, filename=filename)


    