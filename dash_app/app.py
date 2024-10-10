import dash as ds
from Materials_Data_Analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import base64
import io
import pickle


app = ds.Dash(__name__, use_pages=True)
app.config.suppress_callback_exceptions = True
app.title = 'Cyclic Voltammogram Analysis'

app.layout = ds.html.Div([
    ds.html.H1('Materials characterization analytics contents'),
    ds.html.Div([
        ds.html.Div(
            ds.dcc.Link(f"{page['name']} Analysis", href=page["relative_path"])
        ) for page in ds.page_registry.values()
    ]),
    ds.page_container
])

if __name__ == '__main__':
    app.run(debug=True)

def run():
    app.run_server(debug=True)
