import plotly.graph_objects as go
import plotly.express as px

axis_dict = {
    'color': '#85929E',
    'dividercolor': '#625D5D',
    'gridcolor': 'black',
    'linecolor': 'black',
    'showline': False,
    'showgrid': False,
    'linewidth': 0,
    'zeroline': False
}

custom_dark_template = dict(
    layout=go.Layout(
        plot_bgcolor='black',
        paper_bgcolor='black',
        xaxis=axis_dict,
        yaxis=axis_dict,
        margin={'t': 60, 'b': 60, 'l': 60, 'r': 60},
        width=1000,
        height=600,
        colorway=px.colors.qualitative.Set1,
        font={'color': '#EAECEE'}
    )
)
