# <a name="SetupXChemCARforDevelopers"></a>**Setup XChem-CAR for Developers**

Instructions for downloading and running XChem-CAR for developers,<br>

**If you wish to use XChem-CAR you are advised to use the webapp:** <br>
[URL ONCE LIVE]

**Continue with this guide if you wish to download as setup XChem-CAR
for development purposes**<br><br>

## <a name="Setting up Windows Subsystem for Linux"></a>Windows subsystem for linux

If you're using Windows, to install git-crypt, it's strongly advised that you install Windows Subsystem for Linux (WSL).
For seting up WSL2 - you can follow these instructions: https://www.digitalocean.com/community/tutorials/how-to-install-the-windows-subsystem-for-linux-2-on-microsoft-windows-10 and/or https://docs.microsoft.com/en-gb/windows/wsl/install-win10

## <a name="VisualStudioCode"></a>Visual Studio Code

these instructions are designed for Visual Studio Code which can be installed for free from: https://code.visualstudio.com/

## <a name="GitCryptKey"></a>Git-Crypt Key

Secrets required for running CAR are encrypted, to unencrypt and run you will need the key from the XChem-CAR software maintainer<br><br>

# <a name="RepositoryfromGitHub"></a>Clone the "xchem-car" repository from GitHub
In you home directory eg `/home/<username>`, clone  xchem-car repo using:

`git clone https://github.com/Waztom/xchem-CAR.git` <br>

XChem-CAR uses GitHub for version control<br>
to get started with working on CAR clone the <em>"xchem-CAR"</em> Repository from github to your device.<br><br>

### <a name="UsefulGitHubbranches"></a>Useful GitHub branches - pass for setting up

| Branch  | Description                                                | URL                                              |
| ------- | ---------------------------------------------------------- | ------------------------------------------------ |
| Main    | Most recent, stable, release                               | https://github.com/Waztom/xchem-CAR              |
| Develop | new features will be added <br> here before being released | https://github.com/Waztom/xchem-CAR/tree/Develop |

all branches can be found: https://github.com/Waztom/xchem-CAR/branches
<br><br>

# <a name="Docker"></a>Docker

## <a name="InstallDocker"></a>Install Docker

First you'll need Docker Desktop (or the relevent Docker Engine on Linux) you can find the appropriate download like this on https://www.docker.com/get-started

## <a name="InstallDockerCompose"></a>Install Docker Compose

if you're using a Linux machine, once Docker is installed also install docker compose, instructions are available at: https://docs.docker.com/compose/install/

for WSL (Windows), you do not need to install Docker compose

## <a name="InstallVSCodeExtention"></a>Install VS Code Extensions

Docker and Docker Compose should now be installed <br>
<em>(If on Windows/Mac, start docker desktop)</em><br>

Open Visual Studio Code<br>
to check docker is running correctly open the terminal and run:<br>

> `docker --version`<br>

you should get a response similar to:

> Docker version 18.09.2, build 6247962

In Visual Studio Code open the extensions panel (left-hand panel or using Ctrl+Shift+X ) and search for "<em>Remote - Containers extension</em>" and click **Install**.

Once installed a box with two arrows pointing in opposite directions should appear in the bottom left corner of Visual studio code
<br>
<br>

## <a name="InstallRemoteWSL"></a>Install RemoteWSL Extension

Ctrl + SHIFT + X and type in 'Remote - WSL' to install the extension.
See (https://code.visualstudio.com/blogs/2020/07/01/containers-wsl) for more information about using dev containers in WSL2

# <a name="gitcrypt"></a>git-crypt

git-crypt (https://github.com/AGWA/git-crypt) is used for encrypting secrets required to run CAR
you need the appropriate <em>crypt-key</em> file from the software maintainer.

If you are using a Windows machine then it is necessary to download the Windows Subsystem for Linux 2. A very good guide is found here: https://www.digitalocean.com/community/tutorials/how-to-install-the-windows-subsystem-for-linux-2-on-microsoft-windows-10

You will need to allow WSL integration with Docker Desktop. To do this go to settings on Docker > Resources > WSL Integration. Then enable WSL integration for your desired distribution.

Incorporation of VSCode with WSL for further information: https://code.visualstudio.com/docs/remote/wsl-tutorial

**Note: We have tested this using Ubuntu 18.04 and 20.04, compatibility of other Linux distributions have not been investigated.**

## <a name="InstallGitCrypt"></a>Install Git-Crypt

download the compressed git-crypt package (https://www.agwa.name/projects/git-crypt/downloads/git-crypt-0.6.0.tar.gz)

Extract the package to create the directory **git-crypt-0.6.0**:
> `tar -xvzf git-crypt-0.6.0.tar.gz` <br>

in terminal:

> `$ cd git-crypt-0.6.0`<br> >`$ make`<br> >`# make install`

## <a name="UnlockingSecrets"></a>Unlocking Secrets

once Git-Crypt is installed unlock the secrets using:

> `cd [file path to XChem-CARS on your device]`<br> >`git-crypt unlock [path to git-crypt crypt-key]`

# <a name="Startsystem"></a>Start system

### <a name="LocateRepository"></a>Locate Repository

- in terminal change directory to your copy of the repo :

  > `cd `<em>`[local file path to xchem-CAR repository]`

    </em>
    or open VS Code and go to File-> Open Folder and open the repository directory<br>

### <a name="StartRemoteContainer"></a>Start Remote Container

- start Visual Studio remote container with **Ctrl + Shift + P** and type **"Remote-containers: Open folder in container"** then click on that option. <br> ensure you have the repository folder [your file path/xchem-CAR] selected and choose **"Ok"**/**"Open"**
- Your container should start to build, click on the popup notification at the bottom right of visual studio to view the log/progress

### <a name="TimetoLaunch"></a>Time to Launch

- Open a new terminal that you can interact with. if the terminal is visible at the bottom of the screen click on the plus "create new integrated terminal" or use the keybord shortcut "**Ctrl+Shift+`**" button or use the adjacent "split terminal" (or "**Ctrl+Shift+5**") button to see the new terminal adjacent to the current terminal
- you should now be in the container running Debian Linux
- in the new terminal type:
  > `cd CAR` <br> >`mkdir log && touch logsfile.log` <br> >`python3 manage.py makemigrations` <br> >`python3 manage.py migrate` <br> >`npm install --quiet --legacy-peer-deps`<br> >`python3 manage.py runserver`<br>
- to compile the main.js file, in a separate terminal inside your development container:
  > `npm run dev`<br>

If you are only interested in running the application or developing the backend code, you will only need to run the `npm run dev` command
once. For frontend developers, the npm command above tracks any changes made to the frontend code and recompiles the main.js file.

for future launches, you will not need to perform the migrations, install the node packages, compile the main.js and only need to launch the Django server by running: >`cd CAR` <br> >`python3 manage.py runserver`<br>

to upload files in CAR, you need to start a Celery worker in a separate terminal inside your development container:

- open a new terminal the same way as last time ([see Time to Launch](#TimeToLaunch))
- in the new terminal type:
  > `cd CAR`<br> >`celery -A CAR worker -l info`<br>

if you make any changes to the Django models, you will need to run the the migrations again in the CAR directory:

- makemigrations and migrate the Django models:
  > `python3 manage.py makemigrations` <br> >`python3 manage.py migrate` <br>

# Opening the application

at the end of the step "[Time to Launch](#TimeToLaunch)" an address to use the visual interface should have been displayed ("http://127.0.0.1:8000/"), Ctrl+Click on the link in terminal or copy and paste the link into your web browser to use the CAR interface

> http://127.0.0.1:8000/

### <a name="Troubleshooting"></a>Troubleshooting

WLS2/Windows users have reported problems with instructions in "Time to Launch":

### <a name="File Permissions"></a>File Permissions

Sometimes, when a new branch is created on Github, try to `git pull`, get a keyerror for accessing the environment variables in the development container or trying to run the `npm run dev` yielding a `Error: EACCES: permission denied`, you need to change the file permissions on the repo folder.

In your WSL2 Linux terminal **NB outside** your dev container:

- change the file permissions of the repository folder
  `sudo chown -R <username> <path to xchem-CAR>`
