param([int]$Port = 8000)
$env:PYTHONUTF8 = "1"
python "$PSScriptRoot/launch.py" backend --port $Port
