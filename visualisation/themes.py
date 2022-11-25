import plotly.graph_objects as go

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
        colorway=['#E1E1E1'],
        font={'color': '#EAECEE'}
    )
)
