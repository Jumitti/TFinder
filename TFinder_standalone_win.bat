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

:: Check if Python >= 3.10 is installed
echo %COLOR_CYAN%Checking for Python version higher or equal to 3.10...%COLOR_RESET%
for /f "usebackq tokens=2 delims= " %%v in (`python --version`) do set PYTHON_VERSION=%%v

:: Extract major and minor version numbers
for /f "tokens=1,2 delims=." %%a in ("%PYTHON_VERSION%") do (
    set MAJOR=%%a
    set MINOR=%%b
)

:: Check if version is >= 3.10
if not defined MAJOR (
    echo %COLOR_RED%Python is not installed.%COLOR_RESET%
    goto install_python
)

if %MAJOR% lss 3 (
    echo %COLOR_RED%Python version is less than 3.10.%COLOR_RESET%
    goto install_python
) else if %MAJOR% equ 3 if %MINOR% lss 10 (
    echo %COLOR_RED%Python version is less than 3.10.%COLOR_RESET%
    goto install_python
)

echo %COLOR_GREEN%Python version %PYTHON_VERSION% is sufficient.%COLOR_RESET%
goto continue

:install_python
echo %COLOR_YELLOW%Installing Python 3.11.9 from https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe...%COLOR_RESET%
powershell -Command "Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe -OutFile \"python-3.11.9-amd64.exe\""
if exist "python-3.11.9-amd64.exe" (
    echo %COLOR_CYAN%Running Python installer...%COLOR_RESET%
    start /wait "python-3.11.9-amd64.exe" /quiet InstallAllUsers=1 PrependPath=1
    del "python-3.11.9-amd64.exe"
) else (
    echo %COLOR_RED%Failed to download Python installer. Check your internet connection.%COLOR_RESET%
    pause
    exit /b
)

:: Verify the new installation
python --version | findstr /r "3\.[1-9][0-9]*" >nul 2>&1
if %errorlevel% neq 0 (
    echo %COLOR_RED%Python installation failed. Please install it manually by downloading on https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe.%COLOR_RESET%
    pause
    exit /b
)
echo %COLOR_GREEN%Python 3.11.9 installed successfully! %COLOR_RESET%

:continue
:: Virtual environment creation
if not exist ".venv" (
    echo %COLOR_CYAN%Creating Python 3.11.9 virtual environment...%COLOR_RESET%
    python -m venv .venv
    if %errorlevel% neq 0 (
        echo %COLOR_RED%Failed to create the virtual environment. Ensure Python 3.11.9 is installed.%COLOR_RESET%
        pause
        exit /b
    )
)
echo %COLOR_GREEN%Virtual environment created successfully! %COLOR_RESET%

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
pip install -r requirements.txt --verbose
if %errorlevel% neq 0 (
    echo %COLOR_RED%Failed to install required packages. Check the requirements.txt file.%COLOR_RESET%
    pause
    exit /b
)
echo %COLOR_GREEN%Python packages installed successfully! %COLOR_RESET%

:: Run Streamlit/TFinder
echo %COLOR_CYAN%Launching the TFinder Streamlit app... (CTRL + C to shutdown) %COLOR_RESET%
streamlit run TFinder-v1.py
if %errorlevel% neq 0 (
    echo %COLOR_RED%Failed to launch the Streamlit app. Check your script for errors.%COLOR_RESET%
    pause
    exit /b
)

:: Disabling venv on shutdown
deactivate
pause