import dash as ds
from analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import base64
import io

source_options = [{'label': 'Biologic', 'value': 'biologic'}, {'label': 'Aftermath', 'value': 'aftermath'}]
upload_button_style = {'width': '50%', 'height': '60px', 'lineHeight': '60px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'textAlign': 'center'}
table_styles = {'width': '50%', 'overflowX': 'auto'}
plotly_template = 'simple_white'

app = ds.Dash(__name__)
app.config.suppress_callback_exceptions = True

app.layout = [

    ### Select the source of the data ###
    ds.html.H2('Select the source of the data'),
    ds.dcc.Dropdown(id='data_source', options=source_options, placeholder='Select the source of the data', style={'width': '60%'}),
    ds.html.Div(style={'margin-bottom': '10px'}),
    ds.dcc.Input(id='scan_rate_input', type='text', placeholder='Enter the scan rate in mV/s if using Aftermath', style={'width': '35%'}),
    ds.html.Div(style={'margin-bottom': '100px'}),

    ### Upload the data file ###
    ds.html.H2('Upload a data file to see the Cyclic Voltammetry analysis'),
    ds.dcc.Upload(id='data_upload', children=ds.html.Div(['Drag and Drop or ', ds.html.A('Select Files')]), style=upload_button_style),
    ds.html.Div(id='file_name', style={'margin-top': '10px'}),
    ds.html.Div(style={'margin-bottom': '100px'}),

    ### Display the basic analysis ###
    ds.html.H2('Basic Analysis'),
    ds.html.Div(id='raw_data'),
    ds.html.Div(id='current_potential_plot'),
    ds.html.Div(id='current_time_plot'),
    ds.html.Div(id='potential_time_plot'),
    ds.html.Div(style={'margin-bottom': '100px'}),

    ### Display the charged passed analysis ###
    ds.html.H2('Charges Passed'),
    ds.html.Div(id='charge_passed_table'),
    ds.html.Div(style={'margin-bottom': '10px'}),
    ds.html.Div(id='charge_passed_plot')

    ]
    
### update the file name ###
@app.callback(
        ds.Output('file_name', 'children'),
        [ds.Input('data_upload', 'filename')]
)
def update_file_name(filename):
    if filename is not None:
        return f"Uploaded file: {filename}"
    else:
        return "No file uploaded"

### Listen for the data upload for a CV ###
@app.callback(
        [ds.Output('raw_data', 'children'), 
         ds.Output('current_potential_plot', 'children'),
         ds.Output('current_time_plot', 'children'),
         ds.Output('potential_time_plot', 'children'),
         ds.Output('charge_passed_table', 'children'),
         ds.Output('charge_passed_plot', 'children')], 
        [ds.Input('data_upload', 'contents'),
         ds.Input('data_source', 'value'), 
         ds.Input('scan_rate_input', 'value')]
)
def display_cv_data(file_contents, source, scan_rate):

    try: 

        the_cv = CyclicVoltammogram.from_html_base64(file_contents = file_contents, source=source, scan_rate = scan_rate)

        data = the_cv.data.round(7).to_dict('records')
        potential_current_plot = the_cv.get_current_potential_plot(width=1400, height=800, template=plotly_template)
        current_time_plot = the_cv.get_current_time_plot(width=1400, height=600, template=plotly_template)
        potential_time_plot = the_cv.get_potential_time_plot(width=1400, height=600, template=plotly_template)
        charge_passed_data = the_cv.get_charge_passed().round(7).to_dict('records')
        charge_passed_plot = the_cv.get_charge_passed_plot(width=1400, height=600, template=plotly_template)

        data_table_element = ds.dash_table.DataTable(data=data, page_size=10, style_table=table_styles)
        potential_vs_current_plot_element = ds.dcc.Graph(figure=potential_current_plot)
        current_vs_time_plot_element = ds.dcc.Graph(figure=current_time_plot)
        potential_vs_time_plot_element = ds.dcc.Graph(figure=potential_time_plot)
        charge_passed_table_element = ds.dash_table.DataTable(data=charge_passed_data, page_size=10, style_table=table_styles)
        charge_passed_plot_element = ds.dcc.Graph(figure=charge_passed_plot)

        return (data_table_element, potential_vs_current_plot_element, current_vs_time_plot_element, 
                potential_vs_time_plot_element, charge_passed_table_element, charge_passed_plot_element)
    
    except Exception as e:
        return [ds.html.Div(["An error occurred: ", str(e)]) for i in range(6)]
        # return [ds.html.Div() for i in range(6)]

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8051, debug=False)
    # app.run_server(debug=True)
