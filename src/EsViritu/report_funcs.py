from dash import Dash, html
import dash_ag_grid as dag
from datetime import datetime
import json
import os
import polars as pl
import tempfile
import numpy as np

## make json, load json, make report with sparklines
def make_dash_aggrid_html_report(windows: 'pl.DataFrame', meta: 'pl.DataFrame', output_html: str = 'report.html') -> str:
    """
    Create an HTML file with Dash from two polars DataFrames, rendering sparklines using the Sparks-Bar-Medium font.
    - windows: DataFrame with columns including Accession and average_coverage
    - meta: DataFrame with columns Accession, mean_coverage, read_count, family, genus, species, RPKMF
    - output_html: Path to save the HTML file
    
    Returns:
        Path to the generated HTML file
    """
    import logging
    logger = logging.getLogger(__name__)
    
    # Helper: map normalized values to Sparks-Bar-Medium glyphs (8 bars: U+2581 to U+2588)
    def sparkline(vals):
        if not vals or len(vals) == 0:
            return ''
            
        try:
            # Unicode block chars: 0x2581 (lowest) to 0x2588 (highest)
            bars = [chr(code) for code in range(0x2581, 0x2589)]
            arr = np.array(vals, dtype=float)
            
            # Handle NaN values
            arr = np.nan_to_num(arr)
            
            if np.all(arr == arr[0]) or np.max(arr) == np.min(arr):  # flat or all same value
                # Use middle bar for flat data
                return bars[3] * len(arr)
            else:
                # Normalize and map to bars
                arr_norm = (arr - arr.min()) / (arr.max() - arr.min())
                idxs = (arr_norm * (len(bars) - 1)).round().astype(int)
                return ''.join([bars[i] for i in idxs])
        except Exception as e:
            logger.warning(f"Error generating sparkline: {e}")
            return '—' * 10  # Fallback if error occurs

    try:
        # Group windows by Accession, flatten average_coverage to cov_list
        grouped = windows.group_by('Accession').agg(
                pl.col('average_coverage').flatten().alias('cov_list')
            )
        # Join with meta on Accession
        joined = grouped.join(meta, on='Accession', how='inner')
        
        # Select and order columns
        columns = [
            'Accession', 'mean_coverage', 'read_count', 'family', 'genus', 'species', 'RPKMF', 'cov_list'
        ]
        # Only include columns that exist in the joined dataframe
        available_columns = [col for col in columns if col in joined.columns]
        joined = joined.select(available_columns)
        
        # Save to temp JSON file
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.json', delete=False) as tmp:
            json_path = tmp.name
            joined.write_json(json_path)
            
        # Load JSON
        with open(json_path) as f:
            data = json.load(f)
            
        # Update columns to match what's actually available
        columns = available_columns
        
    except Exception as e:
        logger.error(f"Error processing data: {e}")
        raise

    # Build Dash app with a pure HTML table
    from dash import Dash, html, dcc
    app = Dash(__name__, assets_folder='assets')
    
    # Define table styles
    table_style = {
        "width": "100%",
        "borderCollapse": "collapse",
        "fontFamily": "Arial, sans-serif",
        "marginTop": "20px",
        "boxShadow": "0 4px 8px 0 rgba(0,0,0,0.2)"
    }
    
    header_style = {
        "backgroundColor": "#f2f2f2",
        "color": "#333",
        "fontWeight": "bold",
        "textAlign": "left",
        "padding": "12px 15px",
        "borderBottom": "1px solid #ddd"
    }
    
    cell_style = {
        "padding": "10px 15px",
        "borderBottom": "1px solid #ddd"
    }
    
    row_style = {"borderBottom": "1px solid #ddd"}
    alt_row_style = {"backgroundColor": "#f9f9f9", "borderBottom": "1px solid #ddd"}
    
    # Format column headers with more readable names
    column_display_names = {
        "Accession": "Accession",
        "mean_coverage": "Mean Coverage",
        "read_count": "Read Count",
        "family": "Family",
        "genus": "Genus",
        "species": "Species",
        "RPKMF": "RPKM-F",
        "cov_list": "Coverage Profile"
    }
    
    # Table header
    table_header = html.Thead([
        html.Tr([
            html.Th(column_display_names.get(col, col), style=header_style) for col in columns
        ])
    ])
    
    # Table rows with alternating colors
    table_rows = []
    for i, row in enumerate(data):
        row_cells = []
        current_style = alt_row_style if i % 2 else row_style
        
        for col in columns:
            val = row.get(col, "")
            cell_attrs = {"style": {**cell_style}}
            
            # Format numeric values
            if col in ["mean_coverage", "RPKMF"] and val:
                try:
                    val = f"{float(val):.2f}"
                except (ValueError, TypeError):
                    pass
            elif col == "read_count" and val:
                try:
                    val = f"{int(val):,}"
                except (ValueError, TypeError):
                    pass
            
            # Special handling for sparklines
            if col == "cov_list":
                spark = sparkline(val)
                cell = html.Td(
                    spark,
                    style={
                        **cell_style,
                        "fontFamily": "Sparks-Bar-Medium", 
                        "fontSize": "18px", 
                        "letterSpacing": "1px",
                        "backgroundColor": "#fff"
                    },
                    className="bar-medium"
                )
            else:
                cell = html.Td(val, style={**cell_style})
                
            row_cells.append(cell)
        table_rows.append(html.Tr(row_cells, style=current_style))
    table_body = html.Tbody(table_rows)

    # Create a more complete layout with title, timestamp, and styling
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    
    app.layout = html.Div([
        # Header section
        html.Div([
            html.H1("Virus Coverage Report", style={"color": "#2c3e50", "marginBottom": "5px"}),
            html.P(f"Generated: {timestamp}", style={"color": "#7f8c8d", "marginTop": "0"}),
        ], style={"textAlign": "center", "marginBottom": "20px"}),
        
        # CSS link for sparkline fonts
        html.Link(rel="stylesheet", href="/assets/sparks-fonts.css"),
        
        # Main table with improved styling
        html.Table([
            table_header,
            table_body
        ], style=table_style),
        
        # Footer with information
        html.Div([
            html.Hr(style={"margin": "30px 0 20px 0"}),
            html.P([
                "Report generated by EsViritu. ",
                html.Span("Coverage profiles ", style={"fontWeight": "bold"}),
                "show normalized coverage across the genome."
            ], style={"color": "#7f8c8d", "fontSize": "0.9em", "textAlign": "center"})
        ])
    ], style={
        "maxWidth": "1200px", 
        "margin": "0 auto", 
        "padding": "20px",
        "fontFamily": "Arial, sans-serif"
    })

    # Save as HTML (without launching server)
    try:
        # Use a more modern approach to generate standalone HTML
        app.layout = app.layout  # Ensure layout is set
        
        # Method 1: Use dash.dash.get_relative_path if available
        try:
            html_string = app.index_string
            resources = app._collect_and_register_resources()
            css = '\n'.join([
                f'<link rel="stylesheet" href="{resource.relative_package_path}">' 
                for resource in resources.get('css', [])
            ])
            js = '\n'.join([
                f'<script src="{resource.relative_package_path}"></script>' 
                for resource in resources.get('js', [])
            ])
            
            # Add custom CSS for sparkline fonts
            css += '\n<style>\n@font-face {\n  font-family: "Sparks-Bar-Medium";\n  src: url("https://cdn.jsdelivr.net/gh/aftertheflood/sparks@9.0.0/output/assets/fonts/Sparks-Bar-Medium.woff2") format("woff2");\n  font-weight: normal;\n  font-style: normal;\n}\n</style>'
            
            html_string = html_string.replace('</head>', f'{css}</head>')
            html_string = html_string.replace('</body>', f'{js}</body>')
            rendered = html_string.replace('{%app_entry%}', app._generate_scripts_html())
            
        # Method 2: Fallback to simpler approach
        except Exception as inner_e:
            logger.warning(f"Using fallback HTML generation method: {inner_e}")
            from dash import html as dash_html
            
            # Create a basic HTML structure
            rendered = f'''
            <!DOCTYPE html>
            <html>
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1">
                <title>Virus Coverage Report</title>
                <style>
                    @font-face {{
                        font-family: "Sparks-Bar-Medium";
                        src: url("https://cdn.jsdelivr.net/gh/aftertheflood/sparks@9.0.0/output/assets/fonts/Sparks-Bar-Medium.woff2") format("woff2");
                        font-weight: normal;
                        font-style: normal;
                    }}
                    body {{
                        font-family: Arial, sans-serif;
                        max-width: 1200px;
                        margin: 0 auto;
                        padding: 20px;
                    }}
                    table {{
                        width: 100%;
                        border-collapse: collapse;
                        margin-top: 20px;
                        box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2);
                    }}
                    th {{
                        background-color: #f2f2f2;
                        color: #333;
                        font-weight: bold;
                        text-align: left;
                        padding: 12px 15px;
                        border-bottom: 1px solid #ddd;
                    }}
                    td {{
                        padding: 10px 15px;
                        border-bottom: 1px solid #ddd;
                    }}
                    tr:nth-child(even) {{
                        background-color: #f9f9f9;
                    }}
                    .bar-medium {{
                        font-family: "Sparks-Bar-Medium";
                        font-size: 18px;
                        letter-spacing: 1px;
                        background-color: #fff;
                    }}
                </style>
            </head>
            <body>
                {app.layout.to_plotly_json()}
            </body>
            </html>
            '''
        
        # Write the HTML file
        with open(output_html, 'w') as f:
            f.write(rendered)
        
        logger.info(f"HTML report saved to {output_html}")
        
        # Clean up temp files
        try:
            if os.path.exists(json_path):
                os.remove(json_path)
        except Exception as e:
            logger.warning(f"Failed to remove temporary file: {e}")
            
        return output_html
        
    except Exception as e:
        logger.error(f"Failed to generate HTML report: {e}")
        raise


