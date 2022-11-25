import plotly.graph_objects as go

custom_dark_template = dict(
    layout=go.Layout(
        plot_bgcolor='black',
        paper_bgcolor='black',
        xaxis={'color': 'white', 'dividercolor': 'white', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        yaxis={'color': 'white', 'dividercolor': 'white', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        margin={'t': 40, 'b': 40, 'l': 40, 'r': 40},
        width=1000,
        height=600,
        colorway=['#E5E4E2']
    )
)
