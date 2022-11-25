import plotly.graph_objects as go

custom_dark_template = dict(
    layout=go.Layout(
        plot_bgcolor='black',
        paper_bgcolor='black',
        xaxis={'color': '#625D5D', 'dividercolor': '#625D5D', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        yaxis={'color': '#625D5D', 'dividercolor': '#625D5D', 'gridcolor': 'black', 'linecolor': 'black', 'showline': False, 'showgrid': False},
        margin={'t': 40, 'b': 40, 'l': 40, 'r': 40},
        width=1000,
        height=600,
        colorway=['#E5E4E2'],
        font={'color': '#E5E4E2'}
    )
)
