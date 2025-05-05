import argparse
import os
from pathlib import Path
import shutil
import re
import zipfile

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate HTML report')
    parser.add_argument('--sample-paths', nargs='+', required=True, help='Paths to sample directories')
    parser.add_argument('--output-dir', default='html_reports', help='Output directory for the reports')
    return parser.parse_args()

def extract_group_sample_info(path):
    """Extract group and sample ID from path name"""
    path_obj = Path(path)
    dir_name = path_obj.name
    
    # Expecting format like "GROUP-SAMPLEID"
    if '-' in dir_name:
        group_name, sample_id = dir_name.rsplit('-', 1)
        return group_name, sample_id, path_obj
    else:
        # Fallback if format is different
        return dir_name, dir_name, path_obj

def generate_sample_html(sample_id, sample_dir, relative_path):
    """Generate HTML content for a single sample"""
    # Find all report files
    fastqc_initial_files = list(sample_dir.glob('*_initial_*_fastqc.html'))
    fastqc_processed_files = list(sample_dir.glob('*_processed_*_fastqc.html'))
    fastqc_files = fastqc_initial_files + fastqc_processed_files
    fastqc_names = [f"{relative_path}/{f.name}" for f in fastqc_files]
    
    # Find other report files
    fastp_file = next(sample_dir.glob('*.fastp.html'), None) or next(sample_dir.glob('*.adapt.html'), None)
    debridged_file = next(sample_dir.glob('*.codes.tsv_SE.png'), None)
    restr_sites_file = next(sample_dir.glob('*_last_oligos.tsv.Nucl_distr.png'), None)
    
    # Find additional analysis files
    filtered_data_files = list(sample_dir.glob('filtered_data_*.png'))
    filtered_out_files = list(sample_dir.glob('filtered_out_data_*.png'))
    hisat2_logs = list(sample_dir.glob('*.hisat2.summary.log'))
    
    # Create the FastQC buttons HTML
    fastqc_initial_buttons = ""
    fastqc_processed_buttons = ""
    
    for f in fastqc_initial_files:
        display_name = f.name.replace("_fastqc.html", "")
        fastqc_initial_buttons += f'<button class="nav-button" onclick="loadReport(\'{relative_path}/{f.name}\')">{display_name}</button> '
    
    for f in fastqc_processed_files:
        display_name = f.name.replace("_fastqc.html", "")
        fastqc_processed_buttons += f'<button class="nav-button" onclick="loadReport(\'{relative_path}/{f.name}\')">{display_name}</button> '
    
    # Create the navigation buttons HTML
    nav_buttons = ""
    if fastqc_files:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'fastqc-section-' + sample_id + '\')">FastQC Reports</button>'
    if fastp_file:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'fastp-section-' + sample_id + '\')">FastP Report</button>'
    if debridged_file:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'debridged-section-' + sample_id + '\')">Debridged Analysis</button>'
    if restr_sites_file:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'restr-sites-section-' + sample_id + '\')">Restriction Sites</button>'
    if filtered_data_files:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'filtered-data-section-' + sample_id + '\')">Filtered Data Plots</button>'
    if filtered_out_files:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'filtered-out-section-' + sample_id + '\')">Filtered Out Data</button>'
    if hisat2_logs:
        nav_buttons += '<button class="nav-button" onclick="showSection(\'alignment-section-' + sample_id + '\')">Alignment Logs</button>'

    # Generate filtered data image HTML
    filtered_data_html = ""
    for img in filtered_data_files:
        filtered_data_html += f'''
        <div class="image-container">
            <h3>{img.stem.replace("_", " ").title()}</h3>
            <img src="{relative_path}/{img.name}" alt="{img.stem}">
        </div>
        '''
    
    # Generate filtered out data image HTML
    filtered_out_html = ""
    for img in filtered_out_files:
        filtered_out_html += f'''
        <div class="image-container">
            <h3>{img.stem.replace("_", " ").title()}</h3>
            <img src="{relative_path}/{img.name}" alt="{img.stem}">
        </div>
        '''
    
    # Generate HISAT2 logs HTML
    hisat2_html = ""
    for log_file in hisat2_logs:
        hisat2_html += f'<h3>{log_file.name}</h3><pre class="log-content">'
        try:
            with open(log_file, 'r') as f:
                hisat2_html += f.read()
        except:
            hisat2_html += f"Error reading {log_file.name}"
        hisat2_html += '</pre>'

    sample_html = f'''
    <div class="sample-container">
        <h2 id="sample-{sample_id}">{sample_id}</h2>
        
        <div class="nav-buttons">
            {nav_buttons}
        </div>

        <div id="fastqc-section-{sample_id}" class="section">
            <h3>FastQC Reports</h3>
            
            <div class="tab-content">
                <h4>Initial FastQC</h4>
                <div class="nav-buttons">
                    {fastqc_initial_buttons}
                </div>
                
                <h4>Processed FastQC</h4>
                <div class="nav-buttons">
                    {fastqc_processed_buttons}
                </div>
                
                <iframe id="fastqcFrame-{sample_id}" src="{fastqc_names[0] if fastqc_names else ''}"></iframe>
            </div>
        </div>

        <div id="fastp-section-{sample_id}" class="section">
            <h3>FastP Report</h3>
            <iframe src="{relative_path}/{fastp_file.name if fastp_file else ''}" style="width:100%; height:800px;"></iframe>
        </div>

        <div id="debridged-section-{sample_id}" class="section">
            <h3>Debridged Analysis</h3>
            <div class="image-container">
                <img src="{relative_path}/{debridged_file.name if debridged_file else ''}" alt="Debridged analysis">
            </div>
        </div>

        <div id="restr-sites-section-{sample_id}" class="section">
            <h3>Restriction Sites Distribution</h3>
            <div class="image-container">
                <img src="{relative_path}/{restr_sites_file.name if restr_sites_file else ''}" alt="Restriction sites distribution">
            </div>
        </div>
        
        <div id="filtered-data-section-{sample_id}" class="section">
            <h3>Filtered Data Plots</h3>
            {filtered_data_html}
        </div>
        
        <div id="filtered-out-section-{sample_id}" class="section">
            <h3>Filtered Out Data</h3>
            {filtered_out_html}
        </div>
        
        <div id="alignment-section-{sample_id}" class="section">
            <h3>Alignment Statistics</h3>
            {hisat2_html}
        </div>
    </div>
    '''
    return sample_html

def generate_group_html_report(group_name, samples_info, output_dir):
    """Generate a combined HTML report for a group of samples"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create sample navigation
    sample_nav = ""
    for sample_id, _, _ in samples_info:
        sample_nav += f'<button class="group-nav-button" onclick="showSample(\'{sample_id}\')">{sample_id}</button>'
    
    # Generate HTML for each sample
    sample_html_sections = ""
    for sample_id, _, sample_path in samples_info:
        # Create a relative path from the output directory to the sample directory
        rel_path = os.path.basename(sample_path)
        sample_html_sections += generate_sample_html(sample_id, sample_path, rel_path)
    
    html_content = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{group_name} - Group Analysis Report</title>
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
            h1, h2, h3 {{
                color: #333;
            }}
            h1 {{
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
            .group-nav-button {{
                padding: 10px 20px;
                background-color: #2196F3;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                transition: background-color 0.3s;
            }}
            .group-nav-button:hover {{
                background-color: #0b7dda;
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
                padding: 10px;
                border: 1px solid #ddd;
                border-radius: 4px;
            }}
            img {{
                max-width: 100%;
                height: auto;
            }}
            .section {{
                display: none;
                margin-top: 20px;
            }}
            .section.active {{
                display: block;
            }}
            .sample-container {{
                display: none;
                margin-top: 30px;
                border-top: 2px solid #ddd;
                padding-top: 20px;
            }}
            .sample-container.active {{
                display: block;
            }}
            .log-content {{
                background-color: #f8f8f8;
                padding: 15px;
                border-radius: 4px;
                overflow-x: auto;
                white-space: pre-wrap;
                font-size: 14px;
            }}
            .tab-content {{
                margin-top: 15px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>{group_name} - Group Analysis Report</h1>
            
            <h2>Samples</h2>
            <div class="nav-buttons">
                {sample_nav}
            </div>

            {sample_html_sections}
        </div>

        <script>
            // Keep track of currently active section type
            let currentSectionType = '';
            
            function loadReport(reportPath) {{
                // Get the current sample ID from the iframe ID
                const iframeId = document.activeElement.closest('.sample-container').querySelector('iframe').id;
                document.getElementById(iframeId).src = reportPath;
            }}

            function showSection(sectionId) {{
                // Get the parent sample container
                const sampleContainer = document.getElementById(sectionId).closest('.sample-container');
                
                // Hide all sections in this sample container
                sampleContainer.querySelectorAll('.section').forEach(section => {{
                    section.classList.remove('active');
                }});
                
                // Show selected section
                document.getElementById(sectionId).classList.add('active');
                
                // Store the section type (strip sample ID suffix)
                currentSectionType = sectionId.substring(0, sectionId.lastIndexOf('-'));
            }}

            function showSample(sampleId) {{
                // Hide all sample containers
                document.querySelectorAll('.sample-container').forEach(container => {{
                    container.classList.remove('active');
                }});
                
                // Find the sample container by the ID in the h2 tag
                const sampleContainers = document.querySelectorAll('.sample-container');
                let foundContainer = null;
                
                for (let container of sampleContainers) {{
                    if (container.querySelector('h2').id === 'sample-' + sampleId) {{
                        container.classList.add('active');
                        foundContainer = container;
                        break;
                    }}
                }}
                
                if (foundContainer) {{
                    // If we have an active section type and it exists in this sample, show it
                    if (currentSectionType) {{
                        const sectionId = currentSectionType + '-' + sampleId;
                        const section = document.getElementById(sectionId);
                        
                        if (section) {{
                            // Hide all sections in the container
                            foundContainer.querySelectorAll('.section').forEach(s => {{
                                s.classList.remove('active');
                            }});
                            
                            // Show the matching section
                            section.classList.add('active');
                            return;
                        }}
                    }}
                    
                    // Fallback: show the first section if we couldn't find the matching section
                    const firstSection = foundContainer.querySelector('.section');
                    if (firstSection) {{
                        firstSection.classList.add('active');
                    }}
                }}
            }}

            // Initialize the page - show the first sample and its first section
            document.addEventListener('DOMContentLoaded', function() {{
                const firstSampleContainer = document.querySelector('.sample-container');
                if (firstSampleContainer) {{
                    firstSampleContainer.classList.add('active');
                    
                    // Show first section of first sample
                    const firstSection = firstSampleContainer.querySelector('.section');
                    if (firstSection) {{
                        firstSection.classList.add('active');
                        // Initialize the currentSectionType
                        currentSectionType = firstSection.id.substring(0, firstSection.id.lastIndexOf('-'));
                    }}
                }}
            }});
        </script>
    </body>
    </html>
    '''

    output_file = output_path / f"report_{group_name}.html"
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    return output_file

def create_zip_archive(group_name, samples_info, output_dir, html_file):
    """Create a ZIP archive containing the HTML report and all referenced files"""
    output_path = Path(output_dir)
    zip_path = output_path / f"{group_name}.zip"
    
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        # Add the HTML report
        zipf.write(html_file, os.path.basename(html_file))
        
        # Add sample directories with their content
        for _, _, sample_path in samples_info:
            sample_name = os.path.basename(sample_path)
            
            for root, _, files in os.walk(sample_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.join(sample_name, os.path.relpath(file_path, sample_path))
                    zipf.write(file_path, arcname)
    
    return zip_path

def main():
    args = parse_arguments()
    
    # Group samples by group name
    groups = {}
    for sample_path in args.sample_paths:
        group_name, sample_id, path = extract_group_sample_info(sample_path)
        
        if group_name not in groups:
            groups[group_name] = []
            
        groups[group_name].append((sample_id, group_name, path))
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate reports for each group
    for group_name, samples_info in groups.items():
        html_file = generate_group_html_report(group_name, samples_info, output_dir)
        zip_file = create_zip_archive(group_name, samples_info, output_dir, html_file)
        print(f"Created report for {group_name}: {zip_file}")

if __name__ == '__main__':
    main()