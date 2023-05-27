# Modelleren 2B
Google Drive: https://drive.google.com/drive/folders/1t3hshQnNEx3ZVykXBI5LH9K99vD6owvA?usp=share_link

Planning: https://docs.google.com/document/d/183vVBEheg9Vm7GuTYzWlVo_zXa3OJrRZxNFNfCXTA-w/edit?usp=sharing

Overleaf: https://www.overleaf.com/project/644902bedd93e9eb3f3369a1

# Python environement setup
Please ensure you have Python version 3.10 installed.

## Using a normal Python installation
### Using cmd
- Open `cmd.exe` in the root of `modelling-2b` as administrator and run `python -m venv venv && "venv/Scripts/activate.bat" && pip install -r requirements.txt`.
- To activate virtual environment use `"venv/Scripts/activate.bat"`.
- To deactivate virtual environment use `deactivate`.

### Using Git Bash
- Open Git Bash in the root of `modelling-2b` as administrator and run `python -m venv venv && source venv/Scripts/activate && pip install -r requirements.txt`.
- To activate virtual environment use `source venv/Scripts/activate`.
- To deactivate virtual environment use `deactivate`.

### Selecting Python interpreter in Visual Studio Code
- To select a the virtual environment, use the **Python: Select Interpreter** command from the **Command Palette** (<kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>P</kbd>) and enter `.\venv\Scripts\python.exe`.

## Using an Anaconda installation
Replace `C:/ProgramData/Anaconda3/` with `your/path/to/Anaconda3` in the code below in case they do not match.
- Open `cmd.exe` in the root of `modelling-2b` as administrator and run `"C:/ProgramData/Anaconda3/condabin/activate.bat" && conda env create -f environment.yml && conda activate modelling-2b-venv`.
- To activate virtual environment use `conda activate modelling-2b-venv`.
- To deactivate virtual environment use `conda deactivate`.

### Selecting Python interpreter in Anaconda
The interpreter is located at `C:\ProgramData\Anaconda3\envs\modelling-2b-venv\python.exe`.
- To use Spyder, either activate the virtual environment and run `conda install spyder` and run `spyder` or activate the virtual environment from the Spyder console.
- To use Jupyter Notebook, activate the virutal environment and run `jupyter notebook`. Once you have selected a file you can select the virtual environment from **Kernel** > **Change kernel** > **Python 3 (ipykernel)**.
