import dash as ds
from analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import base64
import io
import pickle

source_options = [{'label': 'Biologic', 'value': 'biologic'}, {'label': 'Aftermath', 'value': 'aftermath'}]
upload_button_style = {'width': '50%', 'height': '30px', 'lineHeight': '30px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'textAlign': 'center'}
table_styles = {'width': '50%', 'overflowX': 'auto'}
button_style = {'width': '30%', 'height': '30px', 'lineHeight': '15px', 'borderWidth': '2px', 'borderStyle': 'dashed', 'textAlign': 'center'}
plotly_template = 'ggplot2'

ds.register_page(__name__)

layout = ds.html.Div([

    ds.html.Br(),
    ds.html.Hr(),
    ds.html.Br(),
    ds.html.H1("Welcome to Cyclic Voltammogram Analysis", style={'textAlign': 'center'}),

    ### Input box ###
    ds.html.Div(children=[
        ds.html.Hr(),
        ds.html.H2('Load Cyclic Voltammogram'),
        ds.dcc.Dropdown(id='data_source', options=source_options, placeholder='Select the source of the data', style={'width': '60%'}),
        ds.dcc.Input(id='scan_rate_input', type='text', placeholder='Enter the scan rate in mV/s if using Aftermath', style={'width': '35%'}),
        ds.html.Br(),
        ds.html.Br(),
        ds.html.Br(),
        ds.dcc.Upload(id='data_upload', children=ds.html.Div(['Drag and Drop or Click here to ', ds.html.A('Select File')]), style=upload_button_style),
        ds.html.Div(id='file_name', style={'margin-top': '10px'}),
        ds.html.Br(),
        ds.html.Br(),
        ds.html.Button('Load Cyclic Voltammogram', id='store_cv_button', style=button_style),
        ds.html.Div(id='file_load_status', style={'margin-top': '10px'}),
        ds.dcc.Store(id="cv_stored")
        ], style={'padding': 10, 'flex': 10}),

    ds.html.Hr(),

    ### Display the basic analysis ###
    ds.html.Div(children=[
        ds.html.H2('Basic Analysis'),
        ds.html.Div(id='raw_data'),
        ds.html.Br(),
        ds.html.Div(id='current_potential_plot'),
        ds.html.Div(id='current_time_plot'),
        ds.html.Div(id='potential_time_plot')
        ], style={'padding': 10, 'flex': 10}),

    ds.html.Hr(),

    ### Display the charge passed analysis ###
    ds.html.H2('Charge Passed Analysis'),
    ds.html.Div(id='charge_passed_table'),
    ds.html.Br(),
    ds.html.Div(children=ds.html.Div(id='charge_passed_plot'), style={'display': 'inline-block', 'width': '49%'}),
    ds.html.Div(children=ds.html.Div(id='current_integration_plot_blank'), style={'display': 'inline-block', 'width': '49%'}),

    ds.html.Hr()

    ])
    


def pickle_and_encode(data):
    pickled_data = pickle.dumps(data)
    encoded_data = base64.b64encode(pickled_data).decode('utf-8')
    return encoded_data



def pickle_and_decode(encoded_data):
    pickled_data = base64.b64decode(encoded_data)
    data = pickle.loads(pickled_data)
    return data



@ds.callback(
        ds.Output('file_name', 'children'),
        [ds.Input('data_upload', 'filename')]
)
def update_file_name(filename):
    """
    Callback to update the file selected text
    """
    if filename is not None:
        return f"Selected file: {filename}"
    else:
        return "No file selected"
    


@ds.callback(
    [ds.Output('cv_stored', 'data'),
     ds.Output('file_load_status', 'children')],
    ds.Input('store_cv_button', 'n_clicks'),
    [ds.State('data_upload', 'contents'),
     ds.State('data_source', 'value'), 
     ds.State('scan_rate_input', 'value'),
     ds.State('data_upload', 'filename')],
    prevent_initial_call=True
)
def store_cv_data(n_clicks, file_contents, source, scan_rate, file_name):
    """
    Callback to store the CV data and update the text to let the user know they updated the CV text
    """
    the_cv = CyclicVoltammogram.from_html_base64(file_contents = file_contents, source=source, scan_rate = scan_rate)
    encoded_cv = pickle_and_encode(the_cv)
    return encoded_cv, f"Created Cyclic Voltammogram based off of the {source} file, {file_name} :)"



@ds.callback(
    [ds.Output('raw_data', 'children'),
     ds.Output('current_potential_plot', 'children'),
     ds.Output('current_time_plot', 'children'),
     ds.Output('potential_time_plot', 'children')],
    [ds.Input('cv_stored', 'data')],
    prevent_initial_call=True
)
def display_basic_analysis(encoded_cv):
    """
    Callback to display the basic analysis of the CV
    """
    the_cv = pickle_and_decode(encoded_cv)

    data = the_cv.data.round(7).to_dict('records')
    data_table_element = ds.dash_table.DataTable(data=data, page_size=10, style_table=table_styles)

    current_time_plot = the_cv.get_current_time_plot(width=1400, height=600, template=plotly_template)
    potential_time_plot = the_cv.get_potential_time_plot(width=1400, height=500, template=plotly_template)
    potential_current_plot = the_cv.get_current_potential_plot(width=1100, height=800, template=plotly_template)

    potential_vs_current_plot_element = ds.dcc.Graph(figure=potential_current_plot)
    current_vs_time_plot_element = ds.dcc.Graph(figure=current_time_plot)
    potential_vs_time_plot_element = ds.dcc.Graph(figure=potential_time_plot)

    return data_table_element, potential_vs_current_plot_element, current_vs_time_plot_element, potential_vs_time_plot_element



@ds.callback(
    [ds.Output('charge_passed_table', 'children'),
     ds.Output('charge_passed_plot', 'children'),
     ds.Output('current_integration_plot_blank', 'children')],
    [ds.Input('cv_stored', 'data')],
    prevent_initial_call=True
)
def display_charge_passed_analysis(encoded_cv):
    """
    Callback to display the charge passed analysis
    """
    the_cv = pickle_and_decode(encoded_cv)

    charge_passed_table = the_cv.get_charge_passed().round(7).to_dict('records')
    charge_passed_table_element = ds.dash_table.DataTable(data=charge_passed_table, page_size=10, style_table=table_styles)

    charge_passed_plot = the_cv.get_charge_passed_plot(width=740, height=600, template=plotly_template)
    charge_passed_plot_element = ds.dcc.Graph(figure=charge_passed_plot, id='charge_passed_id')

    interactive_integration_element = ds.dcc.Graph(id='current_integration_plot')

    return charge_passed_table_element, charge_passed_plot_element, interactive_integration_element



@ds.callback(
    ds.Output('current_integration_plot', 'figure'),
    [ds.Input('charge_passed_id', 'clickData')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def update_integration_plot(clickData, encoded_cv):
    """
    Callback to update the integration plot based off of the click data
    """
    the_cv = pickle_and_decode(encoded_cv)

    if clickData is None:
        return the_cv.get_current_time_plot(width=700, height=600, template=plotly_template)

    point_info = clickData['points'][0]
    cycle = point_info['x']
    direction = point_info['customdata'][0]

    integration_plot = the_cv.get_charge_integration_plot(cycle=cycle, direction=direction, width=700, height=600, template=plotly_template)

    return integration_plot
