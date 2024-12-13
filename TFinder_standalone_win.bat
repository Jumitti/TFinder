@echo off
:: Colors
for /f "delims=" %%i in ('echo prompt $E^|cmd') do set "ESC=%%i"
set "COLOR_RESET=%ESC%[0m"
set "COLOR_GREEN=%ESC%[32m"
set "COLOR_RED=%ESC%[31m"
set "COLOR_YELLOW=%ESC%[33m"
set "COLOR_BLUE=%ESC%[34m"
set "COLOR_CYAN=%ESC%[36m"
set "COLOR_BOLD=%ESC%[1m"

:: Welcome
echo %COLOR_GREEN%======================================================================================================%COLOR_RESET%
echo %COLOR_GREEN%!%COLOR_RESET%                                        %COLOR_BOLD%%COLOR_BLUE%Welcome to TFinder%COLOR_RESET%                                         %COLOR_GREEN%!%COLOR_RESET%
echo %COLOR_GREEN%======================================================================================================%COLOR_RESET%
echo %COLOR_CYAN%TFinder%COLOR_RESET% is a %COLOR_BOLD%%COLOR_GREEN%Python easy-to-use web tool%COLOR_RESET% for identifying Transcription Factor Binding Sites (TFBS) and
echo Individual Motif (IM). Using the %COLOR_BOLD%NCBI API%COLOR_RESET%, it can easily extract either the promoter or terminal regions of a gene 
echo through a simple query of NCBI gene name or ID. It enables simultaneous analysis across five different species 
echo for an unlimited number of genes. The tool searches for TFBS and IM in different formats, including IUPAC codes
echo and JASPAR entries. Moreover, %COLOR_CYAN%TFinder%COLOR_RESET% also allows the generation and use of a %COLOR_BOLD%Position Weight Matrix (PWM)%COLOR_RESET%. 
echo Finally, the data may be recovered in a tabular form and a graph showing the relevance of the TFBSs and IMs as 
echo well as its location relative to the Transcription Start Site (TSS) or gene end. The results may be sent by email 
echo to the user facilitating the ulterior analysis and data sharing.
echo %COLOR_GREEN%======================================================================================================%COLOR_RESET%
echo %COLOR_YELLOW%Created by Minniti Julien%COLOR_RESET% - %COLOR_BLUE%GitHub(%COLOR_RESET%https://github.com/Jumitti/TFinder%COLOR_BLUE%)%COLOR_RESET%
echo %COLOR_YELLOW%MIT Licence%COLOR_RESET% - %COLOR_BLUE%https://github.com/Jumitti/TFinder/blob/main/LICENSE%COLOR_RESET%
echo %COLOR_GREEN%======================================================================================================%COLOR_RESET%

:: Check if Python 3.11.9 is installed
echo %COLOR_CYAN%Checking for Python 3.11.9 installation...%COLOR_RESET%
python --version | find "3.11.9" >nul 2>&1
if %errorlevel% neq 0 (
    echo %COLOR_YELLOW%Python 3.11.9 is not installed.%COLOR_RESET%
    echo %COLOR_CYAN%Downloading and installing Python 3.11.9...%COLOR_RESET%

    :: Download the Python 3.11.9 installer (Windows x64)
    powershell -Command "Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe -OutFile python-3.11.9-amd64.exe"
    if exist python-3.11.9-amd64.exe (
        echo %COLOR_CYAN%Running Python installer...(may be long)%COLOR_RESET%
        start /wait python-3.11.9-amd64.exe /quiet InstallAllUsers=1 PrependPath=1
        del python-3.11.9-amd64.exe
    ) else (
        echo %COLOR_RED%Failed to download Python installer. Check your internet connection.%COLOR_RESET%
        pause
        exit /b
    )
    
    :: Check again if Python 3.11.9 is installed
    python --version | find "3.11.9" >nul 2>&1
    if %errorlevel% neq 0 (
        echo %COLOR_RED%Python 3.11.9 installation failed. Please install it manually.%COLOR_RESET%
        pause
        exit /b
    )
    echo %COLOR_GREEN%Python 3.11.9 installed successfully! %COLOR_RESET%
)

:: Virtual environment creation
if not exist ".venv" (
    echo %COLOR_CYAN%Creating Python 3.11.9 virtual environment...%COLOR_RESET%
    python -m venv .venv
    if %errorlevel% neq 0 (
        echo %COLOR_RED%Failed to create the virtual environment. Ensure Python 3.11.9 is installed.%COLOR_RESET%
        pause
        exit /b
    )
    echo %COLOR_GREEN%Virtual environment created successfully! %COLOR_RESET%
)

:: Activation venv
echo %COLOR_CYAN%Activating virtual environment...%COLOR_RESET%
call .venv\Scripts\activate
if %errorlevel% neq 0 (
    echo %COLOR_RED%Failed to activate the virtual environment.%COLOR_RESET%
    pause
    exit /b
)
echo %COLOR_GREEN%Virtual environment activated successfully! %COLOR_RESET%

:: Updating pip and installing dependencies
echo %COLOR_CYAN%Installing/updating Python packages...%COLOR_RESET%
pip install --upgrade pip
pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo %COLOR_RED%Failed to install required packages. Check the requirements.txt file.%COLOR_RESET%
    pause
    exit /b
)
echo %COLOR_GREEN%Python packages installed successfully! %COLOR_RESET%

:: Run Streamlit/TFinder
echo %COLOR_CYAN%Launching the TFinder Streamlit app... (CTRL + C to shutdown)%COLOR_RESET%
streamlit run TFinder-v1.py
if %errorlevel% neq 0 (
    echo %COLOR_RED%Failed to launch the Streamlit app. Check your script for errors.%COLOR_RESET%
    pause
    exit /b
)

:: Disabling venv on shutdown
deactivate
pause