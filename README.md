---
title: Library Explorer
emoji: 🦀
colorFrom: pink
colorTo: gray
sdk: gradio
sdk_version: 6.10.0
app_file: app.py
pinned: false
license: mit
---

# ExploDIA
### A Spectral Library Analysis and Visualization Tool for Proteomics

ExploDIA is an interactive tool designed for biologists and proteomics researchers to analyze, merge, filter, and visualize spectral libraries — no coding experience required.

---

## 📋 Table of Contents
- [What is ExploDIA?](#what-is-explodia)
- [Two Ways to Use ExploDIA](#two-ways-to-use-explodia)
- [Getting Started (Downloadable Version)](#getting-started-downloadable-version)
- [Features](#features)
- [File Requirements](#file-requirements)
- [Frequently Asked Questions](#frequently-asked-questions)

---

## What is ExploDIA?

ExploDIA is a web-based tool that runs in your browser. It allows you to:
- View key statistics about your spectral library
- Merge two libraries together
- Filter a library by modifications or proteins
- Visualize precursor ion mobility vs m/z space
- Explore peptides for a specific protein of interest

No installation of Python, R, or any other programming language is required.

---

## Two Ways to Use ExploDIA

### 🌐 Online Version (Hugging Face)
**Link:** [huggingface.co/spaces/kpotula02/library-explorer](https://huggingface.co/spaces/kpotula02/library-explorer)

- Access directly from any browser
- No downloads or setup required
- Best for files **under ~1.5 GB**
- May be slower on very large files due to server memory limits

### 💻 Downloadable Version (Recommended for Large Files)
**Link:** [github.com/kapishp/library-explorer](https://github.com/kapishp/library-explorer)

- Runs entirely on your own computer
- Your data never leaves your machine — ideal for sensitive research data
- Significantly faster file loading and processing
- **Recommended for files above ~1.5 GB**
- Works without an internet connection after setup

---

## Getting Started (Downloadable Version)

### Step 1 — Check that Python is installed
ExploDIA requires Python 3 to be installed on your computer.

**Mac:**
1. Open the **Terminal** app (search for it using Spotlight: press `Cmd + Space` and type "Terminal")
2. Type the following and press Enter:
3. If you see a version number (e.g. `Python 3.11.0`), you are good to go
4. If you see an error, download and install Python from [python.org](https://python.org)

**Windows:**
1. Open the **Command Prompt** (press `Windows key`, type "cmd", press Enter)
2. Type the following and press Enter:
3. If you see a version number, you are good to go
4. If you see an error, download and install Python from [python.org](https://python.org). Make sure to check **"Add Python to PATH"** during installation

---

### Step 2 — Download ExploDIA
1. Go to [github.com/kapishp/library-explorer](https://github.com/kapishp/library-explorer)
2. Find the file called `LibraryExplorer.zip`
3. Click on it and then click **Download**
4. Once downloaded, find the zip file (usually in your Downloads folder) and unzip it
   - **Mac:** Double click the zip file
   - **Windows:** Right click the zip file and select "Extract All"

---

### Step 3 — Launch the App

**Mac:**
1. Open the `LibraryExplorer` folder
2. Double click `run_mac.command`
3. If you see a security warning saying the file "cannot be verified":
   - Go to **System Settings** → **Privacy & Security**
   - Scroll down and click **Open Anyway**
4. A terminal window will open and the app will start automatically in your browser after about 15 seconds
5. If the browser shows an error at first, wait a few more seconds and refresh the page

**Windows:**
1. Open the `LibraryExplorer` folder
2. Double click `run_windows.bat`
3. If Windows shows a security warning, click **More Info** → **Run Anyway**
4. A command prompt window will open and the app will start automatically in your browser after about 15 seconds

> ⚠️ **Important:** Do not close the terminal/command prompt window while using the app — it is what keeps the app running. To stop the app, simply close that window.

> ⏱️ **First launch:** The very first time you run ExploDIA it will automatically install the required packages. This may take 3–5 minutes. After that, it will start up quickly every time.

---

## Features

### 📊 Tab 1 — Library Explorer
Analyze the statistics of your spectral library.

**How to use:**
1. Click **Upload Library File** and select your `.tsv` file
2. Wait for the file to load (this may take a few minutes for large files)
3. Use the **Filter by Protein** dropdown to view stats for a specific protein, or select **Overall** for the full library

**What you will see:**
- Key statistics: total rows, unique proteins, unique peptides, retention time range, ion mobility range
- Interactive charts:
  - Frequency vs Normalized Retention Time
  - Frequency vs Precursor m/z
  - Average Retention Time vs m/z
  - Average Retention Time vs Frequency
  - Precursor Charge vs Precursor m/z (bubble chart)

---

### 🔗 Tab 2 — Library Merger
Merge two spectral libraries into one combined file.

**How to use:**
1. Upload your **First Library** and **Second Library** TSV files
2. Choose how to handle duplicates:
   - **Remove duplicates from Library 1** — any entry in Library 1 that also appears in Library 2 will be removed before merging (recommended)
   - **Keep all duplicates** — all entries from both libraries are kept
3. Choose how to handle column differences between the two files:
   - **Keep only common columns** — only columns present in both files are kept
   - **Preserve all columns** — all columns from both files are kept, with empty values where data is missing
4. Click **Merge Libraries**
5. Once complete, click **Download Merged Library** to save the result

**What you will see:**
- A detailed log showing how many rows each library contained, how many duplicates were found, and the final row count after merging

---

### 🧹 Tab 3 — Library Extractor
Filter your spectral library to create a smaller, customized version.

> ⚠️ **Note:** This tab works best with files under ~500 MB. For larger files, processing may take several minutes.

**How to use:**
1. Upload your library file
2. Wait for the file to load — you will see a list of detected UniMod modifications and proteins
3. **Step 2 — Remove Modifications:** Select any modifications you want to remove from the library and click **Apply Mod Filter**
4. **Step 3 — Filter Proteins:** Choose whether to keep only specific proteins, remove specific proteins, or skip protein filtering. Use **Find Related Isoforms** to detect protein variants
5. **Step 4 — Export:** Enter a filename and click **Finalize & Download** to save your filtered library

---

### 📡 Tab 4 — Method Visualization
Visualize your spectral library in m/z vs Ion Mobility space, broken down by charge state.

**How to use:**
1. Upload your library file
2. Optionally upload a **PASEF Window File** if you want to overlay acquisition windows on the plot
3. If uploading a PASEF file, select the correct **PASEF Method Type** from the dropdown:
   - DIA-PASEF (rectangles)
   - Diagonal-PASEF (polygons)
   - Slice Multi Window
   - Slice Simple (Cycle 1)
4. Choose an **Analysis Mode**:
   - **None (just plot)** — shows all precursors in grey
   - **Map UniMod Modifications** — highlights specific modifications in color. Enter the UniMod IDs you want to highlight (e.g. `21,35`) or type `ALL` to highlight all detected modifications
   - **Map Missed Cleavages** — colors points by the number of missed cleavages
5. Optionally check **Remove duplicate coordinates before plotting** to deduplicate points with the same m/z and ion mobility values
6. Click **Generate Visualization**

**What you will see:**
- Six interactive scatter plots: one showing all precursors, and one for each charge state (+1 through +5)
- Hover over any point to see the peptide sequence, protein, charge, m/z, ion mobility, and PTM information
- PASEF windows overlaid in red (if uploaded)
- You can customize the colors of PTM highlights using the **Customize PTM Colors** panel

> 💡 **Tip:** The plot shows a representative sample of up to 200,000 points for performance. The full dataset is still used for all statistics.

---

### 🔬 Tab 5 — Protein Peptide Cloud
Visualize all peptides for a specific protein of interest on top of the full library.

**How to use:**

**Step 1 — Load Library & Select Protein**
1. Upload your library file
2. Wait for the protein list to populate
3. Select a protein from the **Select Protein ID** dropdown, or type a protein ID directly (e.g. `P10636`)

**Step 2 — PASEF Overlay (Optional)**
1. Upload a PASEF window file if desired
2. Check **Overlay PASEF Windows** and select the correct PASEF method type

**Step 3 — Generate Overview Plot**
1. Click **Generate Peptide Cloud (Plot 1)**
2. The plot will show:
   - All library precursors as small grey dots (background)
   - Your selected protein's peptides highlighted on top:
     - **Green circles** = unmodified peptides
     - **Colored diamonds** = modified peptides (color indicates modification type)
3. Hover over any point to see full peptide details
4. Download the full precursor table for your protein using the **Download** button

**Step 4 — Highlight Specific Peptides**
1. Use the **Select Peptide(s) to Highlight** dropdown to choose one or more peptides from the list
   - 🔶 indicates a modified peptide
   - ⚪ indicates an unmodified peptide
2. Alternatively, enter serial numbers manually in the **Sl_No** field (e.g. `1,5,12`)

**Step 5 — Generate Highlighted Plot**
1. Click **Generate Highlighted Plot (Plot 2)**
2. The selected peptides will be labeled and highlighted with colored markers
3. A table showing the details of the selected peptides will appear below the plot

---

## File Requirements

| Tab | Required File | Format |
|-----|--------------|--------|
| Library Explorer | Spectral library | `.tsv`, `.csv`, or `.txt` |
| Library Merger | Two spectral libraries | `.tsv`, `.csv`, or `.txt` |
| Library Extractor | Spectral library | `.tsv`, `.csv`, or `.txt` |
| Method Visualization | Spectral library | `.tsv`, `.csv`, or `.txt` |
| Method Visualization | PASEF window file (optional) | `.txt` or `.csv` |
| Protein Peptide Cloud | Spectral library | `.tsv`, `.csv`, or `.txt` |
| Protein Peptide Cloud | PASEF window file (optional) | `.txt` or `.csv` |

**Required columns in spectral library files:**
- `PrecursorMz`
- `PrecursorCharge`
- `PrecursorIonMobility`
- `PeptideSequence`
- `ModifiedPeptideSequence`
- `ProteinId`
- `NormalizedRetentionTime`

---

## Frequently Asked Questions

**How long will it take to load my file?**
This depends on your file size and whether you are using the online or downloadable version. As a rough guide:
- Under 500 MB: 1–2 minutes
- 500 MB – 1.5 GB: 2–5 minutes
- Above 1.5 GB: Use the downloadable version for best results

**My file is above 1.5 GB — will it work?**
Yes, but we strongly recommend using the downloadable version rather than the online version. The downloadable version uses your computer's own memory and is significantly faster for large files.

**The app opened in my browser but shows an error — what do I do?**
Wait 10–15 seconds and refresh the page. The app may still be starting up. If the error persists, close the terminal window and launch the app again.

**Can I use ExploDIA offline?**
Yes! The downloadable version works completely offline once it has been set up for the first time.

**Is my data safe?**
- **Online version:** Files are temporarily uploaded to Hugging Face's servers for processing and are not stored permanently
- **Downloadable version:** Your data never leaves your computer at any point

**The Library Extractor (Tab 3) is taking a very long time — is that normal?**
Tab 3 needs to scan every row of your file to find all unique proteins and modifications, which takes longer than other tabs. For files above 500 MB this may take several minutes. We are working on improving this in a future update.

**Can I use any PASEF window file?**
Yes — the PASEF window file does not need to correspond to the specific library file. It represents your instrument's acquisition settings and can be reused across different libraries acquired with the same method.

**What UniMod IDs should I enter for common modifications?**

| Modification | UniMod ID |
|-------------|-----------|
| Phosphorylation | 21 |
| Oxidation | 35 |
| Carbamidomethyl | 4 |
| Acetylation | 1 |
| Deamidation | 7 |
| Ubiquitinyl | 121 |

---

## Support
For questions or issues, please open a GitHub Issue at [github.com/kapishp/library-explorer/issues](https://github.com/kapishp/library-explorer/issues)
