process GenerateReport {

    // Input txt file containing statistics
    input:
    tuple val(meta), path(stats_file)

    // Output html report file
    output:
    file '*.html', emit: html_report

    script:
    """
    # Parse the statistics from the input file and generate the HTML report with separate tabs for each sample

    python3 <<EOF
    from jinja2 import Template
    import json

    with open('${stats_file}', 'r') as f:
        stats = json.load(f)

    template = Template('<html><body>{% for sample, data in stats.items() %}<h1>{{ sample }}</h1><table>{% for key, value in data.items() %}<tr><td>{{ key }}</td><td>{{ value }}</td></tr>{% endfor %}</table>{% endfor %}</body></html>')
    rendered = template.render(stats=stats)

    with open('report_${meta.prefix}.html', 'w') as f:
        f.write(rendered)
    EOF
    """
}