#!/bin/bash
cd "$(dirname "$0")"

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    osascript -e 'display alert "Python 3 is not installed. Please install it from python.org first." as critical'
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Setting up Library Explorer for the first time..."
    python3 -m venv venv
    source venv/bin/activate
    pip install gradio pandas datasets matplotlib numpy plotly
else
    source venv/bin/activate
fi

# Run the app
echo "Starting Library Explorer..."
python3 app.py &

# Wait for app to start then open browser
sleep 15

open "http://127.0.0.1:7860"

wait