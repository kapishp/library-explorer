import subprocess
import sys
import os
import webbrowser
import time
import threading

def open_browser():
    time.sleep(8)
    webbrowser.open("http://127.0.0.1:7860")

if __name__ == "__main__":
    threading.Thread(target=open_browser).start()
    os.system(f"{sys.executable} app.py")