import argparse
import os
from pathlib import Path
import shutil

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate HTML report')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    parser.add_argument('--out-dir', required=True, help='Output directory')
    parser.add_argument('--output', required=True, help='Output HTML file')
    return parser.parse_args()

def generate_html_report(sample_id, out_dir, output_file):
    # Get all files in the output directory
    out_dir_path = Path(out_dir)
    
    # Find FastQC files
    fastqc_files = list(out_dir_path.glob('*_fastqc.html'))
    fastqc_names = [f.name for f in fastqc_files]
    
    # Find other files
    fastp_file = next(out_dir_path.glob('*.fastp.html'), None)
    debridged_file = next(out_dir_path.glob('*.codes.tsv_SE.png'), None)
    restr_sites_file = next(out_dir_path.glob('*_last_oligos.tsv.Nucl_distr.png'), None)

    html_content = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{sample_id} - Analysis Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f0f0f0;
            }}
            .container {{
                max-width: 1200px;
                margin: 20px auto;
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 0 10px rgba(0,0,0,0.1);
            }}
            h1, h2 {{
                color: #333;
                text-align: center;
            }}
            .nav-buttons {{
                display: flex;
                flex-wrap: wrap;
                gap: 10px;
                margin-bottom: 20px;
                justify-content: center;
            }}
            .nav-button {{
                padding: 10px 20px;
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                transition: background-color 0.3s;
            }}
            .nav-button:hover {{
                background-color: #45a049;
            }}
            iframe {{
                width: 100%;
                height: 800px;
                border: none;
                border-radius: 4px;
            }}
            .image-container {{
                text-align: center;
                margin: 20px 0;
            }}
            img {{
                max-width: 100%;
                height: auto;
            }}
            .section {{
                display: none;
            }}
            .section.active {{
                display: block;
            }}
        </style>
    </head>
    <body>
        <h1>{sample_id} - Analysis Report</h1>
        
        <div class="nav-buttons">
            <button class="nav-button" onclick="showSection('fastqc-section')">FastQC Reports</button>
            <button class="nav-button" onclick="showSection('fastp-section')">FastP Report</button>
            <button class="nav-button" onclick="showSection('debridged-section')">Debridged Analysis</button>
            <button class="nav-button" onclick="showSection('restr-sites-section')">Restriction Sites</button>
        </div>

        <div id="fastqc-section" class="container section">
            <h2>FastQC Reports</h2>
            <div class="nav-buttons">
                {' '.join([f'<button class="nav-button" onclick="loadReport(\'{name}\')">{name.replace("_fastqc.html", "")}</button>' for name in fastqc_names])}
            </div>
            <iframe id="fastqcFrame" src="{fastqc_names[0] if fastqc_names else ''}"></iframe>
        </div>

        <div id="fastp-section" class="container section">
            <h2>FastP Report</h2>
            <iframe src="{fastp_file.name if fastp_file else ''}" style="width:100%; height:800px;"></iframe>
        </div>

        <div id="debridged-section" class="container section">
            <h2>Debridged Analysis</h2>
            <div class="image-container">
                <img src="{debridged_file.name if debridged_file else ''}" alt="Debridged analysis">
            </div>
        </div>

        <div id="restr-sites-section" class="container section">
            <h2>Restriction Sites Distribution</h2>
            <div class="image-container">
                <img src="{restr_sites_file.name if restr_sites_file else ''}" alt="Restriction sites distribution">
            </div>
        </div>

        <script>
            function loadReport(reportPath) {{
                document.getElementById('fastqcFrame').src = reportPath;
            }}

            function showSection(sectionId) {{
                // Hide all sections
                document.querySelectorAll('.section').forEach(section => {{
                    section.style.display = 'none';
                }});
                // Show selected section
                document.getElementById(sectionId).style.display = 'block';
            }}

            // Show FastQC section by default
            document.getElementById('fastqc-section').style.display = 'block';
        </script>
    </body>
    </html>
    '''

    with open(output_file, 'w') as f:
        f.write(html_content)

def main():
    args = parse_arguments()
    generate_html_report(args.sample_id, args.out_dir, args.output)

if __name__ == '__main__':
    main()