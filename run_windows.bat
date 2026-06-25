@echo off
cd /d "%~dp0"

:: Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo Python is not installed. Please install it from python.org first.
    pause
    exit
)

:: Create virtual environment if it doesn't exist
if not exist "venv" (
    echo Setting up Library Explorer for the first time...
    python -m venv venv
    call venv\Scripts\activate
    pip install gradio pandas datasets matplotlib numpy plotly
) else (
    call venv\Scripts\activate
)

:: Run the app
echo Starting Library Explorer...
start /b python app.py

:: Wait and open browser
timeout /t 6 /nobreak
start http://127.0.0.1:7860