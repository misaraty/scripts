## **[中文版本]()**

## AGPD

**Automatically Generating Python Dependencies (AGPD)** is a single-file Python utility that scans `import` statements in a Python project, identifies the third-party packages actually used by the source code, retrieves their installed versions from the current Python environment, and automatically generates a `requirements.txt` file.

## Usage

AGPD uses Python's built-in `ast` module to parse `.py` files and extract imports such as:

```python
import numpy
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
```

The program automatically:

1. Scans Python files in the project;

2. Extracts imported modules;

3. Excludes Python standard-library modules;

4. Excludes local project modules;

5. Uses `importlib.metadata` to map import names to installed distribution names;

6. Retrieves package versions from the current Python environment;

7. Generates `requirements.txt` and a dependency analysis report.

For example:

```text
cv2      -> opencv-python
PIL      -> Pillow
sklearn  -> scikit-learn
```

The generated dependency file may look like:

```text
numpy==2.2.6
pandas==2.3.1
scikit-learn==1.7.0
```

Package versions are obtained from the Python environment in which AGPD is executed.

## Usage

Place `AGPD.py` in the project directory:

```text
project/
├── AGPD.py
├── main.py
├── model.py
└── utils.py
```

Run it inside the Python, virtualenv, or conda environment used by the project:

```bash
python AGPD.py
```

By default, AGPD scans the current directory and generates:

```text
requirements.txt
requirements_report.txt
```

Scan a specific directory:

```bash
python AGPD.py /path/to/project
```

Scan a single Python file:

```bash
python AGPD.py main.py
```

Specify the requirements output file:

```bash
python AGPD.py . -o requirements.txt
```

Specify the report file:

```bash
python AGPD.py . --report requirements_report.txt
```

After generating `requirements.txt`, install the dependencies in a new environment with:

```bash
pip install -r requirements.txt
```
