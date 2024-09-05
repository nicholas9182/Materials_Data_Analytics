import dash as ds
from analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram

my_cv = CyclicVoltammogram.from_biologic(path = '../test_trajectories/cyclic_voltammetry/biologic5.txt')

app = ds.Dash()

app.layout = [
    ds.html.H2('Select the source of the data'),
    ds.dcc.Dropdown(id='data_source', options=[{'label': 'Biologic', 'value': 'biologic'}, {'label': 'Aftermath', 'value': 'aftermath'}], 
                    value='hide', placeholder='Select the source of the data'),
    ds.html.Div(style={'margin-bottom': '10px'}),
    
    ds.html.Div(id='scan_rate_input_container'),
    ds.html.Div(style={'margin-bottom': '100px'}),


    ds.html.H2('Upload a data file to see the Cyclic Voltammetry analysis'),
    ds.dcc.Upload(id='data_upload', children=ds.html.Div(['Drag and Drop or ', ds.html.A('Select Files')]), multiple=False, 
                  style={'width': '100%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 
                         'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px'}),
    ds.html.Div(style={'margin-bottom': '100px'}),
    
    # ds.dash_table.DataTable(id='raw_data', data=[], columns=[]),
    # ds.dcc.Graph(id='current_potential_plot'),
    # ds.dcc.Graph(id='current_time_plot'),
    # ds.dcc.Graph(id='potential_time_plot'),
    # ds.dash_table.DataTable(id='charge_passed', data=[]),
    # ds.dcc.Graph(id='charge_passed_plot')
    ds.html.H2('Basic Analysis'),
    ds.dash_table.DataTable(data=my_cv.data.round(7).to_dict('records'), page_size=10, style_table={'width': '50%', 'overflowX': 'auto'}),
    ds.dcc.Graph(figure=my_cv.get_current_potential_plot(width=1400, height=800, template='simple_white')),
    ds.dcc.Graph(figure=my_cv.get_current_time_plot(width=1400, height=600, template='simple_white')),
    ds.dcc.Graph(figure=my_cv.get_potential_time_plot(width=1400, height=600, template='simple_white')),
    ds.html.Div(style={'margin-bottom': '100px'}),

    ds.html.H2('Charges Passed'),
    ds.dash_table.DataTable(data=my_cv.get_charge_passed().round(7).to_dict('records'), page_size=10, style_table={'width': '50%', 'overflowX': 'auto'}),
    ds.html.Div(style={'margin-bottom': '60px'}),
    ds.dcc.Graph(figure=my_cv.get_charge_passed_plot(width=1400, height=600, template='simple_white'))
    ]

@app.callback(ds.Output('scan_rate_input_container', 'children'), [ds.Input('data_source', 'value')])
def display_scan_rate_box(selected_value):
    if selected_value == 'aftermath':
        return ds.dcc.Input(id='scan_rate_input', type='text', placeholder='Enter the scan rate in mV/s')
    else:
        return ds.html.Div()

if __name__ == '__main__':
    app.run_server(debug=True)
