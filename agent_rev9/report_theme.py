"""Shared dashboard styling for live and standalone reports."""
import base64
import json
import re
from pathlib import Path


ASSETS = Path(__file__).parent / 'assets'


def style_report(document):
    if 'id="omicraft-report-theme"' in document:
        return document
    logo = base64.b64encode((ASSETS / 'omicraft-logo.png').read_bytes()).decode('ascii')
    brand = (
        '<header class="report-brand" id="omicraft-report-brand">'
        '<a href="/omicraft-version-7-1" aria-label="OmiCraft workspace">'
        '<img src="data:image/png;base64,' + logo + '" alt="OmiCraft by team OMICS2DRUG"></a></header>'
        '<nav class="report-nav" aria-label="Main navigation"><div class="report-nav-inner">'
        '<div class="report-nav-tabs"><a href="/omicraft-version-7-1">Workspace</a>'
        '<a href="#" aria-current="page">Pipeline Report</a></div>'
        '<button type="button" class="report-print" id="reportPrint">Print / PDF</button>'
        '</div></nav>'
    )
    script = (ASSETS / 'report-theme.js').read_text().replace('__REPORT_BRAND__', json.dumps(brand))
    theme = '<style id="omicraft-report-theme">' + (ASSETS / 'report-theme.css').read_text() + '</style>'
    theme += '<script>' + script + '</script>'
    document = re.sub(r'<body\b', '<body data-omicraft-report', document, count=1, flags=re.IGNORECASE)
    return document.replace('</head>', theme + '</head>', 1)
