import dash as ds
from Materials_Data_Analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import base64
import io
import pickle

source_options = [{'label': 'Biologic', 'value': 'biologic'}, {'label': 'Aftermath', 'value': 'aftermath'}]
upload_button_style = {'width': '50%', 'height': '30px', 'lineHeight': '30px', 'borderWidth': '1px', 'borderStyle': 'dashed', 'textAlign': 'center'}
table_styles = {'width': '50%', 'overflowX': 'auto'}
button_style = {'width': '30%', 'height': '30px', 'lineHeight': '15px', 'borderWidth': '2px', 'borderStyle': 'dashed', 'textAlign': 'center'}
slider_style = {'width':'70%'}
plotly_template = 'ggplot2'


ds.register_page(__name__)


introduction_text = """
This page is designed to help you analyse cyclic voltammograms. There are three main concepts here which will help with the analysis:

1. The cycling voltammogram is broken into cycles, sweep direction and segments. 
2. The sweep direction is the direction of the potential sweep, either reduction or oxidation. 
3. The cycles go from a potential peak/trough to the next potential peak/trough in such a way that it will contain exactly one reduction sweep and one oxidation sweep.
4. A segment is an integer denoting a single sweep in a given cycle. Each cycle+direction combination will have a segment 0, 1, 2, etc.


Things to note:

 - When constructing the cyclic voltammogram, current-time roots are added into the data, as these are useful for many types of analysis.
 - The cycles start and end at different points compared to the default for the potentiostat. This is also to aid analysis.
 - This is a tool, and ultimately you are responsible for the analysis and ensuring that the data is correct.

 
Main contributions:
 1. Nicholas Siemons (Chueh Group at Stanford University, Salleo Group at Stanford University).  Email: nsiemons@stanford.edu

 
Contributions:
1. Gerwin Dijk (Salleo Group at Stanford University). Email: gdijk@stanford.edu



Have fun!
"""


layout = ds.html.Div([

    ds.html.Br(), ds.html.Hr(), ds.html.Hr(), ds.html.Br(),
    ds.html.H1("Welcome to cyclic voltammogram analysis", style={'textAlign': 'center'}),
    ds.html.Div(introduction_text, style={'textAlign': 'left', 'whiteSpace': 'pre-line', 'padding': 10}),

    ### Input box ###
    ds.html.Div(children=[
        ds.html.Hr(), ds.html.Hr(),
        ds.html.H2('Load Cyclic Voltammogram'),
        ds.dcc.Dropdown(id='data_source', options=source_options, placeholder='Select the source of the data', style={'width': '60%'}),
        ds.dcc.Input(id='scan_rate_input', type='text', placeholder='Enter the scan rate in mV/s if using Aftermath', style={'width': '35%'}),
        ds.html.Br(), ds.html.Br(), ds.html.Br(),
        ds.dcc.Upload(id='data_upload', children=ds.html.Div(['Drag and Drop or Click here to ', ds.html.A('Select File')]), style=upload_button_style),
        ds.html.Div(id='file_name', style={'margin-top': '10px'}),
        ds.html.Br(), ds.html.Br(),
        ds.html.Button('Get editing parameters for Cyclic Voltammogram', id='get_cv_parameters', style=button_style),
        ds.html.Div(id='parameters_message'),
        ds.dcc.Store(id="cv_parameters_for_editing"),
        ds.dcc.Store(id="cv_stored"),
        ds.html.Br(), ds.html.Br(), ds.html.Br(),
        ds.html.H3('Select the cycles to analyze'),
        ds.html.Div(ds.dcc.RangeSlider(min=0, max=5, value=[0, 5], step=1, marks={i: str(i) for i in range(6)}), id='cycle_slider2', style=slider_style),
        ds.html.Br(),
        ds.html.H3('Down-sample the CV, select the number of potential steps per cycle'),
        ds.html.Br(),
        ds.html.Div(ds.dcc.Slider(min=1, max=1000, value=1000, tooltip={"placement": "bottom", "always_visible": True}), id='downsample_slider2', style=slider_style),
        ds.html.Br(),ds.html.Br(),
        ds.html.Button('Update Cyclic Voltammogram', id='update_cv_button', style=button_style),
        ds.html.Div(id='update_cv_code'),
        ds.html.Br(), ds.html.Br(), 
        ds.html.Div(id='create_cv_code'),
        ds.html.Hr(), ds.html.Hr()
        ], style={'padding': 10, 'flex': 10}),

    ds.dcc.Tabs([
        ds.dcc.Tab(label='Basic Analysis', value='basic_analysis', children=[ds.html.Div(id='basic_analysis')]),
        ds.dcc.Tab(label='Charges passed analysis', value='charge_passed_analysis', children=[ds.html.Div(id='charge_passed_analysis')]),
        ds.dcc.Tab(label='Peak fitting analysis', value='peak_fitting_analysis', children=[ds.html.Div(id='peak_fitting_tools'), ds.html.Div(id='peak_fitting_analysis')]),
    ], value='basic_analysis'),

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
    [ds.Output('cv_parameters_for_editing', 'data'),
     ds.Output('parameters_message', 'children')],
     ds.Input('get_cv_parameters', 'n_clicks'),
    [ds.State('data_upload', 'contents'),
     ds.State('data_source', 'value'), 
     ds.State('scan_rate_input', 'value'),
     ds.State('data_upload', 'filename')],
    prevent_initial_call=True
)
def store_cv_parameteres_for_editing(n_clicks, file_contents, source, scan_rate, file_name):
    """
    Callback to store the CV data and update the text to let the user know they updated the CV text
    """
    the_cv = CyclicVoltammogram.from_html_base64(file_contents = file_contents, source=source, scan_rate = scan_rate)
    max_cycle = the_cv.max_cycle
    potential_steps_per_cycle = the_cv.steps_per_cycle

    cv_data = {'cv': pickle_and_encode(the_cv),
               'file_contents': file_contents,
               'source': source,
               'scan_rate': scan_rate, 
               'max_cycle': max_cycle, 
               'potential_steps_per_cycle': potential_steps_per_cycle
               }

    if source == 'biologic':
        code_snippet = f"""```cyclic_voltammogram = CyclicVoltammogram.from_biologic(file_path='{file_name}')```"""
    else:
        code_snippet = f"""```cyclic_voltammogram = CyclicVoltammogram.from_aftermath(file_path='{file_name}', scan_rate={scan_rate})```"""

    code_snippet_element = ds.dcc.Markdown(code_snippet)

    return cv_data, code_snippet_element


@ds.callback(
    ds.Output('downsample_slider2', 'children'),
    ds.Output('cycle_slider2', 'children'),
    [ds.Input('cv_parameters_for_editing', 'data')],
    prevent_initial_call=True
)
def display_editing_tools(cv_data):
    """
    Callback to display the basic analysis of the CV
    """
    max_cycle = cv_data['max_cycle']
    potential_steps_per_cycle = cv_data['potential_steps_per_cycle']
    max_downsample = potential_steps_per_cycle

    cycle_slider = ds.html.Div(ds.dcc.RangeSlider(id='cycle_slider', min=0, max=max_cycle, value=[0, max_cycle], 
                                                  step=1, marks={i: str(i) for i in range(max_cycle+1)}), style=slider_style)

    downsample_slider = ds.html.Div(ds.dcc.Slider(id='downsample_slider', min=1, max=max_downsample, value=potential_steps_per_cycle, 
                                                  tooltip={"placement": "bottom", "always_visible": True}), style=slider_style)

    return downsample_slider, cycle_slider


@ds.callback(
    [ds.Output('cv_stored', 'data'),
    ds.Output('update_cv_code', 'children')],
    ds.Input('update_cv_button', 'n_clicks'),
    [ds.State('cycle_slider', 'value'),
     ds.State('downsample_slider', 'value'),
     ds.State('cv_parameters_for_editing', 'data')],
    prevent_initial_call=True
)
def update_cv_data(n_clicks, cycle_range, down_sample_n, cv_data):
    """
    Callback to update the CV data based off of the sliders
    """
    the_cv = pickle_and_decode(cv_data['cv'])
    max_cycle = cv_data['max_cycle']
    potential_steps_per_cycle = cv_data['potential_steps_per_cycle']
    
    if cycle_range != [0, max_cycle] and down_sample_n == potential_steps_per_cycle:
        cycles = [i for i in range(cycle_range[0], cycle_range[1]+1)]
        new_cv = the_cv.drop_cycles(keep=cycles)
        code_snippet = f"""```cyclic_voltammogram = cyclic_voltammogram.drop_cycles(keep={cycles})```"""
    elif cycle_range == [0, max_cycle] and down_sample_n != potential_steps_per_cycle:
        new_cv = the_cv.downsample(n=down_sample_n)
        code_snippet = f"""```cyclic_voltammogram = cyclic_voltammogram.downsample(n={down_sample_n})```"""
    elif cycle_range == [0, max_cycle] and down_sample_n == potential_steps_per_cycle:
        new_cv = the_cv
        code_snippet = f"""```cyclic_voltammogram = cyclic_voltammogram```"""
    elif cycle_range != [0, max_cycle] and down_sample_n != potential_steps_per_cycle:
        cycles = [i for i in range(cycle_range[0], cycle_range[1]+1)]
        new_cv = the_cv.drop_cycles(keep=cycles).downsample(n=down_sample_n)
        code_snippet = f"""```cyclic_voltammogram = cyclic_voltammogram.drop_cycles(keep={cycles}).downsample(n={down_sample_n})```"""

    code_snippet_element = ds.dcc.Markdown(code_snippet)

    return pickle_and_encode(new_cv), code_snippet_element


@ds.callback(
    ds.Output('basic_analysis', 'children'),
    [ds.Input('cv_stored', 'data')],
    prevent_initial_call=True
)
def display_basic_analysis(encoded_cv):
    """
    Callback to display the basic analysis of the CV
    """
    the_cv = pickle_and_decode(encoded_cv)

    data = the_cv.data.round(7).to_dict('records')
    current_time_plot = the_cv.get_current_time_plot(width=1400, height=600, template=plotly_template)
    potential_time_plot = the_cv.get_potential_time_plot(width=1400, height=500, template=plotly_template)
    potential_current_plot = the_cv.get_current_potential_plot(width=1100, height=800, template=plotly_template)    

    data_table_element = ds.html.Div([
        ds.html.Br(),
        ds.html.H3('Data post processing'),
        ds.html.Div(["""In the following table, the data is shown post processing. This includes the current, potential, cycle, sweep direction and segment."""]),
        ds.dcc.Markdown(f"""```cyclic_voltammogram.data.round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=data, page_size=10, style_table=table_styles),
        ds.html.Button("Download Table as CSV", id="download_raw_data_button"),
        ds.dcc.Download(id="download_raw_data"),
        ds.html.Br()
        ])

    potential_vs_current_plot_element = ds.html.Div([
        ds.html.Br(),
        ds.html.H3('Current vs Potential plot'),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_current_potential_plot()```"""),
        ds.dcc.Graph(figure=potential_current_plot),
        ds.html.Button("Download figure as PDF", id="download_current_potential_button"),
        ds.dcc.Download(id="download_current_potential"),
        ])

    current_vs_time_plot_element = ds.html.Div([
        ds.html.Br(),
        ds.html.H3('Current vs Time plot'),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_current_time_plot()```"""),
        ds.dcc.Graph(figure=current_time_plot),
        ds.html.Button("Download figure as PDF", id="download_current_time_button"),
        ds.dcc.Download(id="download_current_time")
        ])

    potential_vs_time_plot_element = ds.html.Div([
        ds.html.Br(),
        ds.html.H3('Potential vs Time plot'),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_potential_time_plot()```"""),
        ds.dcc.Graph(figure=potential_time_plot),
        ds.html.Button("Download figure as PDF", id="download_potential_time_button"),
        ds.dcc.Download(id="download_potential_time")
    ])

    basic_analysis_element = ds.html.Div(children=[
        data_table_element,
        ds.html.Br(),
        potential_vs_current_plot_element,
        current_vs_time_plot_element,
        potential_vs_time_plot_element
        ], style={'padding': 10, 'flex': 10})

    return basic_analysis_element


@ds.callback(
    ds.Output('charge_passed_analysis', 'children'),
    [ds.Input('cv_stored', 'data')],
    prevent_initial_call=True
)
def display_charge_passed_analysis(encoded_cv):
    """
    Callback to display the charge passed analysis
    """
    the_cv = pickle_and_decode(encoded_cv)

    charge_passed_table_summary = the_cv.get_charge_passed(average_segments = True).round(7).to_dict('records')
    charge_passed_table = the_cv.get_charge_passed().round(7).to_dict('records')
    charge_passed_plot = the_cv.get_charge_passed_plot(width=740, height=600, template=plotly_template)
    max_charges_passed_table_summary = the_cv.get_maximum_charges_passed(average_sections = True).round(7).to_dict('records')
    max_charges_passed_table = the_cv.get_maximum_charges_passed().round(7).to_dict('records')
    max_charge_passed_plot = the_cv.get_maximum_charge_passed_plot(width=740, height=600, template=plotly_template)

    charge_passed_table_summary_element = ds.html.Div([
        ds.html.H3('Charges passed per cycle summary'),
        ds.html.Div(["""In the following table, the charges passed in the cyclic voltammogram are sectioned by cycle and sweep direction. They are then
                     averaged across the cycles. The error represents the standard error across the cycles."""]),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_charge_passed(average_segments = True).round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=charge_passed_table_summary, style_table=table_styles),
        ds.html.Button("Download Table as CSV", id="charges_summary_download_button"),
        ds.dcc.Download(id="charges_summary_download"),
        ds.html.Br()
        ])

    charge_passed_table_element = ds.html.Div([
        ds.html.Br(),
        ds.html.H3('Charges passed per cycle'),
        ds.html.Div(["""In this table, the data is shown per cycle."""]),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_charge_passed().round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=charge_passed_table, page_size=10, style_table=table_styles),
        ds.html.Button("Download Table as CSV", id="charges_cycle_download_button"),
        ds.dcc.Download(id="charges_cycle_download"),
        ds.html.Br()
        ])

    charge_passed_plot_element = ds.html.Div([
        ds.html.Div([
            ds.dcc.Markdown("""```cyclic_voltammogram.get_charge_passed_plot()```"""),
            ds.dcc.Graph(figure=charge_passed_plot, id='charge_passed_id')
            ]),
        ds.html.Div([
            ds.html.Div(id='current_integration_plot_code'),
            ds.dcc.Graph(id='current_integration_plot')
            ])
        ], style={'display': 'flex', 'justify-content': 'space-around'})
    
    max_charges_passed_table_summary_element = ds.html.Div([
        ds.html.H3('Maximum charges passed summary'),
        ds.html.Div(["""In the following table, the charges passed in the cyclic voltammogram are sectioned by when the current goes from positive to negative or 
                     vice versa. These are called 'sections'. They are then averaged across the cycles. The error represents the standard error across the sections."""]),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_maximum_charges_passed(average_sections = True).round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=max_charges_passed_table_summary, style_table=table_styles),
        ds.html.Button("Download Table as CSV", id="max_charges_summary_download_button"),
        ds.dcc.Download(id="max_charges_summary_download"),
        ds.html.Br(), ds.html.Br()
        ])
    
    max_charges_passed_table_element = ds.html.Div([
        ds.html.H3('Maximum charges passed'),
        ds.html.Div(["""In the following analysis, the current-time curve is divided in sections whereby a new section is started whenever the 
                     current flips from positive to negative or vice versa. By integrating these areas, we get the maximum charges that pass in either
                     the anodic or cathodic direction. The maximum charge passed in each section is shown in the table below."""]),
        ds.dcc.Markdown("""```cyclic_voltammogram.get_maximum_charges_passed().round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=max_charges_passed_table, style_table=table_styles),
        ds.html.Button("Download Table as CSV", id="max_charges_passed_download_button"),
        ds.dcc.Download(id="max_charges_passed_download"),
        ds.html.Br()
        ])
    
    max_charge_passed_plot_element = ds.html.Div([
        ds.html.Div([
            ds.dcc.Markdown("""```cyclic_voltammogram.get_maximum_charge_passed_plot()```"""),
            ds.dcc.Graph(figure=max_charge_passed_plot, id='max_charge_passed_id')
            ]),
        ds.html.Div([
            ds.html.Div(id='max_current_integration_plot_code'),
            ds.dcc.Graph(id='max_current_integration_plot')
            ])
        ], style={'display': 'flex', 'justify-content': 'space-around'})

    charge_analysis_element = ds.html.Div(children=[
        charge_passed_table_summary_element,
        charge_passed_table_element,
        ds.html.Br(), ds.html.Br(),
        ds.html.H3('Charges passed per cycle plot'),
        charge_passed_plot_element,
        ds.html.Br(),
        max_charges_passed_table_summary_element,
        max_charges_passed_table_element,
        ds.html.Br(), ds.html.Br(),
        ds.html.H3('Charges passed per cycle plot'),
        max_charge_passed_plot_element
        ], style={'padding': 2, 'flex': 10})

    return charge_analysis_element


@ds.callback(
    ds.Output('peak_fitting_tools', 'children'),
    [ds.Input('cv_stored', 'data')],
    prevent_initial_call=True
)
def display_peak_fitting_tools(encoded_cv):
    """
    Callback to display the peak fitting analysis
    """
    peak_fitting_tools_element = ds.html.Div([
        ds.html.H3('Peak fitting tools'),
        ds.html.Div(["""In the following section, you can fit a polynomial to the data around a peak. The polynomial is then used to find the peak x and y value."""]),
        ds.html.Br(),
        ds.html.H3("Select the order of the polynomial to fit around the peak"),
        ds.html.Div(ds.dcc.Slider(min=2, max=8, step=2, value=4, id='polynomial_order_slider'), style=slider_style),
        ds.html.Br(),
        ds.html.H3("Select the potential window around the peak to fit your polynomial"),
        ds.html.Div(ds.dcc.Slider(min=0.001, max=0.04, value=0.01, id='window_slider'), style=slider_style),
        ds.html.Br(),
        ds.html.Button('Update Peak Fitting', id='update_peak_fitting', style=button_style),
        ])

    return peak_fitting_tools_element


@ds.callback(
    ds.Output('peak_fitting_analysis', 'children'),
    [ds.Input('update_peak_fitting', 'n_clicks')],
    [ds.State('cv_stored', 'data'),
     ds.State('polynomial_order_slider', 'value'),
     ds.State('window_slider', 'value')],
    prevent_initial_call=True
)
def display_peak_fitting_analysis(n_clicks, encoded_cv, polynomial_order, window):

    the_cv = pickle_and_decode(encoded_cv)
    peak_points = the_cv.get_peaks(window=window, polynomial_order=polynomial_order, summary=True).round(7).to_dict('records')
    peak_plot_reduction = the_cv.get_peak_plot(direction='reduction', window = window, polynomial_order = polynomial_order, width=700, height=500)
    peak_plot_oxidation = the_cv.get_peak_plot(direction='oxidation', window = window, polynomial_order = polynomial_order, width=700, height=500)
    current_figure, potential_figure = the_cv.get_plots_peaks_with_cycle(polynomial_order=polynomial_order, window=window, width=700, height=500)

    peak_points_table_summary_element = ds.html.Div([
        ds.html.H3('Peaks for each cycle'),
        ds.html.Div(["""The following table shows the peak potentials and currents for each cycle. The polynomial order and window size are used to fit a polynomial to the data around the peak."""]),
        ds.dcc.Markdown(f"""```cyclic_voltammogram.get_peaks(window={window}, polynomial_order={polynomial_order}, summary=True).round(7)```"""),
        ds.html.Br(),
        ds.dash_table.DataTable(data=peak_points, style_table=table_styles),
        ds.html.Br()
        ])
    
    peak_plot_element = ds.html.Div([
        ds.html.Div([
            ds.dcc.Markdown(f"""```cyclic_voltammogram.get_peak_plot(direction='reduction', window = {window}, polynomial_order = {polynomial_order})```"""),
            ds.dcc.Graph(figure=peak_plot_reduction, id='peak_plot_reduction_plot')
            ]),
        ds.html.Div([
            ds.dcc.Markdown(f"""```cyclic_voltammogram.get_peak_plot(direction='oxidation', window = {window}, polynomial_order = {polynomial_order})```"""),
            ds.dcc.Graph(figure=peak_plot_oxidation, id='peak_plot_oxidation_plot')
            ])
        ], style={'display': 'flex', 'justify-content': 'space-around'})
    
    peak_position_plot_element = ds.html.Div([
        ds.html.Div([
            ds.dcc.Markdown(f"""```cyclic_voltammogram.get_plots_peaks_with_cycle(polynomial_order={polynomial_order}, window={window})[0]```"""),
            ds.dcc.Graph(figure=current_figure, id='current_peak_plot_plot')
            ]),
        ds.html.Div([
            ds.dcc.Markdown(f"""```cyclic_voltammogram.get_plots_peaks_with_cycle(polynomial_order={polynomial_order}, window={window})[1]```"""),
            ds.dcc.Graph(figure=potential_figure, id='potential_peak_oxidation_plot')
            ])
        ], style={'display': 'flex', 'justify-content': 'space-around'})
    
    peaks_analysis_element = ds.html.Div(children=[
        ds.html.Br(),
        peak_points_table_summary_element,
        ds.html.Br(),
        peak_plot_element,
        ds.html.Br(),
        peak_position_plot_element,
        ds.html.Br()
        ], style={'padding': 2, 'flex': 10})
    
    return peaks_analysis_element


@ds.callback(
    [ds.Output('current_integration_plot', 'figure'),
     ds.Output('current_integration_plot_code', 'children')],
    [ds.Input('charge_passed_id', 'clickData')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def update_integration_plot(clickData, encoded_cv):
    """
    Callback to update the integration plot based off of the click data
    """
    the_cv = pickle_and_decode(encoded_cv)

    point_info = clickData['points'][0]
    cycle = point_info['x']
    direction = point_info['customdata'][0]

    integration_plot = the_cv.get_charge_integration_plot(cycle=cycle, direction=direction, width=700, height=600, template=plotly_template)
    code_snippet = f"""```cyclic_voltammogram.get_charge_integration_plot(cycle={cycle}, direction={direction})```"""
    code_snippet_element = ds.dcc.Markdown(code_snippet)

    return integration_plot, code_snippet_element


@ds.callback(
    [ds.Output('max_current_integration_plot', 'figure'),
     ds.Output('max_current_integration_plot_code', 'children')],
    [ds.Input('max_charge_passed_id', 'clickData')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def update_max_charges_integration_plot(clickData, encoded_cv):
    """
    Callback to update the integration plot based off of the click data
    """
    the_cv = pickle_and_decode(encoded_cv)

    point_info = clickData['points'][0]
    section = point_info['x']

    integration_plot = the_cv.get_maximum_charge_integration_plot(section=section, width=700, height=600, template=plotly_template)
    code_snippet = f"""```cyclic_voltammogram.get_maximum_charge_integration_plot(section={section})```"""
    code_snippet_element = ds.dcc.Markdown(code_snippet)

    return integration_plot, code_snippet_element


@ds.callback(
    ds.Output('charges_summary_download', 'data'),
    [ds.Input('charges_summary_download_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_charges_summary(n_clicks, encoded_cv):
    """
    Callback to download the charges summary table
    """
    the_cv = pickle_and_decode(encoded_cv)
    data = the_cv.get_charge_passed(average_segments = True)
    return ds.dcc.send_data_frame(data.to_csv, "charges_summary.csv")


@ds.callback(
    ds.Output('charges_cycle_download', 'data'),
    [ds.Input('charges_cycle_download_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_charges_cycle(n_clicks, encoded_cv):
    """
    Callback to download the charges cycle table
    """
    the_cv = pickle_and_decode(encoded_cv)
    data = the_cv.get_charge_passed()
    return ds.dcc.send_data_frame(data.to_csv, "charges_cycle.csv")


@ds.callback(
    ds.Output('download_raw_data', 'data'),
    [ds.Input('download_raw_data_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_raw_data(n_clicks, encoded_cv):
    """
    Callback to download the raw data
    """
    the_cv = pickle_and_decode(encoded_cv)
    data = the_cv.data
    return ds.dcc.send_data_frame(data.to_csv, "raw_data.csv")


@ds.callback(
    ds.Output('download_current_potential', 'data'),
    [ds.Input('download_current_potential_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_current_potential(n_clicks, encoded_cv):
    """
    Callback to download the current potential plot
    """
    the_cv = pickle_and_decode(encoded_cv)
    fig = the_cv.get_current_potential_plot(width=1100, height=800, template=plotly_template)
    pdf_file = "/tmp/current_potential.pdf"
    fig.write_image(pdf_file, format='pdf')

    return ds.dcc.send_file(pdf_file, "current_potential.pdf")


@ds.callback(
    ds.Output('download_current_time', 'data'),
    [ds.Input('download_current_time_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_current_time(n_clicks, encoded_cv):
    """
    Callback to download the current time plot
    """
    the_cv = pickle_and_decode(encoded_cv)
    fig = the_cv.get_current_time_plot(width=1400, height=600, template=plotly_template)
    pdf_file = "/tmp/current_time.pdf"
    fig.write_image(pdf_file, format='pdf')
    return ds.dcc.send_file(pdf_file, "current_time.pdf")


@ds.callback(
    ds.Output('download_potential_time', 'data'),
    [ds.Input('download_potential_time_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_potential_time(n_clicks, encoded_cv):
    """
    Callback to download the potential time plot
    """
    the_cv = pickle_and_decode(encoded_cv)
    fig = the_cv.get_potential_time_plot(width=1400, height=500, template=plotly_template)
    pdf_file = "/tmp/potential_time.pdf"
    fig.write_image(pdf_file, format='pdf')
    return ds.dcc.send_file(pdf_file, "potential_time.pdf")


@ds.callback(
    ds.Output('max_charges_passed_download', 'data'),
    [ds.Input('max_charges_passed_download_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_max_charges_passed(n_clicks, encoded_cv):
    """
    Callback to download the max charges passed table
    """
    the_cv = pickle_and_decode(encoded_cv)
    data = the_cv.get_maximum_charges_passed()
    return ds.dcc.send_data_frame(data.to_csv, "max_charges_passed.csv")


@ds.callback(
    ds.Output('max_charges_summary_download', 'data'),
    [ds.Input('max_charges_summary_download_button', 'n_clicks')],
    [ds.State('cv_stored', 'data')],
    prevent_initial_call=True
)
def download_max_charges_summary(n_clicks, encoded_cv):
    """
    Callback to download the max charges summary table
    """
    the_cv = pickle_and_decode(encoded_cv)
    data = the_cv.get_maximum_charges_passed(average_sections = True)
    return ds.dcc.send_data_frame(data.to_csv, "max_charges_summary.csv")


