import plotly.graph_objects as go

custom_dark_template = dict(
    layout=go.Layout(
        plot_bgcolor='black',
        paper_bgcolor='black',
        xaxis={'color': '#85929E', 'dividercolor': '#625D5D', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        yaxis={'color': '#85929E', 'dividercolor': '#625D5D', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        margin={'t': 60, 'b': 60, 'l': 60, 'r': 60},
        width=1000,
        height=600,
        colorway=['#E5E4E2'],
        font={'color': '#EAECEE'}
    )
)
