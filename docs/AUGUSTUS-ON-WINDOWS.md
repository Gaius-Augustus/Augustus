# Building and using AUGUSTUS on Windows 10

To use AUGUSTUS on Windows, we recommend the latest version of the Windows Subsystem for Linux (WSL 2). It provides a GNU/Linux environment - including most command-line tools, utilities, and applications - directly on Windows. (See: <https://docs.microsoft.com/en-us/windows/wsl/about>)


## Windows 10 Version

For the WSL 2 a Windows 10, updated to version 2004, Build 19041 or higher, is required. To check the current windows version enter the `ver` command in Windows Command Prompt.

## Windows Features

To use the WSL 2 on your Windows machine the optional Windows features "Windows Subsystem for Linux" and  "Virtual Machine Platform" have to be enabled. To do this open PowerShell as Administrator and run the following commands.

To run PowerShell with administrative privileges click **Start**, type **PowerShell**, right-click **Windows PowerShell**, and then click **Run as administrator**.

    dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
    dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart

Afterwards, restart your machine to complete the WSL installation.

## Installing a Linux distribution

Before the desired Linux distribution can be installed, the update package for the WSL 2 Linux  kernel  must be installed.  The package can be found here: <https://docs.microsoft.com/en-us/windows/wsl/wsl2-kernel>.

For Augustus the WSL version 2 is recommended. Thus, the default WSL version should be set to 2 with the following  PowerShell command 

    wsl --set-default-version 2

If you are using several Linux distributions with both WSL versions in parallel you can set the  WSL version manually for each distribution as described here: <https://docs.microsoft.com/en-us/windows/wsl/install-win10>.

As Linux distribution for Augustus we recommend a current Ubuntu version.  The following link will open the Microsoft store page for Ubuntu 20.04: <https://www.microsoft.com/de-de/p/ubuntu-2004-lts/9n6svws3rx71?rtc=1&activetab=pivot:overviewtab>. Please install this distribution.  After the successful installation a new entry will be created in the start menu. To initialize the new distribution, start it and follow the instructions. You must set a user name and password after the first startup.

Check the installed Linux distributions and their versions with the following  PowerShell command.

    wsl --list --verbose

If Ubuntu is successfully installed it should look like this.

      NAME            STATE           VERSION
    * Ubuntu-20.04    Running         2

Afterwards, you can follow the documentation for build AUGUSTUS [README](../README.md) and run AUGUSTUS [RUNNING-AUGUSTUS](RUNNING-AUGUSTUS.md).

## General information on usage

To view files in a folder you can use the Windows Explorer. To do so, you can navigate to the corresponding folder and then call:

    explorer.exe .

As an editor, [Visual Studio Code](<https://code.visualstudio.com/>) with its extensive extensions is very well suited in combination with the WSL 2. It is recommended to open files or folders from the WSL command prompt. For example:

    code file_to_open.txt.
