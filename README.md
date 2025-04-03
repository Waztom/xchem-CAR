# <a name="LocalDeploymentXChemCAR"></a>**Local deployment of XChem-CAR**
#
Instructions for downloading and running XChem-CAR locally,<br>

**If you wish to use XChem-CAR you are advised to use the webapp:** <br>
[URL ONCE LIVE]

**Continue with this guide if you wish to download as setup XChem-CAR
for running locally or development purposes**<br><br>

## <a name="Setting up Windows Subsystem for Linux"></a>Windows subsystem for linux

If you're using Windows, to install git-crypt, it's strongly advised that you install Windows Subsystem for Linux (WSL).
For setting up WSL2 - you can follow these instructions: https://www.digitalocean.com/community/tutorials/how-to-install-the-windows-subsystem-for-linux-2-on-microsoft-windows-10 and/or https://docs.microsoft.com/en-gb/windows/wsl/install-win10

## <a name="VisualStudioCode"></a>Visual Studio Code

These instructions are designed for Visual Studio Code which can be installed for free from: https://container.visualstudio.com/

## <a name="GitCryptKey"></a>Git-Crypt Key

Secrets required for running CAR are encrypted, to unencrypt and run you will need the key from the XChem-CAR software maintainer<br><br>

# <a name="RepositoryfromGitHub"></a>Clone the "xchem-car" repository from GitHub
If you do not have git installed, inside a terminal:

`sudo apt install git` <br>

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

First you'll need Docker Desktop (or the relevant Docker Engine on Linux) you can find the appropriate download at: https://www.docker.com/get-started and specifically for Ubuntu, instructions at: https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-20-04

## <a name="InstallDockerCompose"></a>Install Docker Compose

if you're using a Linux machine, once Docker is installed also install docker compose, instructions are available at: https://docs.docker.com/compose/install/

for WSL (Windows), you do not need to install Docker compose

## <a name="InstallVSCodeExtention"></a>Install VSCode Extensions

Docker and Docker Compose should now be installed <br>
<em>(If on Windows/Mac, start docker desktop)</em><br>

Open Visual Studio Code<br>
to check docker is running correctly open the terminal and run:<br>

> `docker --version`<br>

you should get a response similar to:

> Docker version 18.09.2, build 6247962

In Visual Studio Code open the extensions panel (left-hand panel or using Ctrl+Shift+X ) and search for "<em>Remote - Containers</em>" and click **Install**.

once installed a box with two arrows pointing in opposite directions should appear in the bottom left corner of Visual studio code
<br>
<br>

## <a name="InstallRemoteWSL"></a>Install Remote - WSL Extension
If you are running WSL - you need to install the Remote - WSL extension. Skip this step if you're running Ubuntu/Linux. 

Ctrl + SHIFT + X and type in 'Remote - WSL' to install the extension.
See (https://container.visualstudio.com/blogs/2020/07/01/containers-wsl) for more information about using dev containers in WSL2

# <a name="gitcrypt"></a>git-crypt

git-crypt (https://github.com/AGWA/git-crypt) is used for encrypting secrets required to run CAR
you need the appropriate <em>crypt-key</em> file from the software maintainer.

If you are using a Windows machine then it is necessary to download the Windows Subsystem for Linux 2. A very good guide is found here: https://www.digitalocean.com/community/tutorials/how-to-install-the-windows-subsystem-for-linux-2-on-microsoft-windows-10

You will need to allow WSL integration with Docker Desktop. To do this go to settings on Docker > Resources > WSL Integration. Then enable WSL integration for your desired distribution.

Incorporation of VSCode with WSL for further information: https://container.visualstudio.com/docs/remote/wsl-tutorial

**Note: We have tested this using Ubuntu 18.04 and 20.04, compatibility of other Linux distributions have not been investigated.**

## <a name="InstallGitCrypt"></a>Install Git-Crypt
if you are using Ubuntu or Debian, you can install git-crypt by:
>`sudo apt-get update` <br>
>`sudo apt-get install git-crypt` <br>

if you are using a Mac, you can install git-crypt using HomeBrew:
>`brew install git-crypt` <br>

## <a name="UnlockingSecrets"></a>Unlocking Secrets

once Git-Crypt is installed unlock the secrets using:

> `cd xchem-car`<br> 
> `git-crypt unlock <'path to git-crypt crypt-key'>`<br>

# <a name="Startsystem"></a>Start system

### <a name="Start VSCode"></a>Start VSCode (WSL)

  > `code .` <br>

### <a name="Start VSCode"></a>Start VSCode (Ubuntu)

  Open VSCode and go to File-> Open Folder and open the repository directory<br>

### <a name="StartRemoteContainer"></a>Start Remote Container

- start Visual Studio remote container with **Ctrl + Shift + P** and type **"Remote-containers: Open folder in container"** then click on that option. <br> ensure you have the repository folder [your file path/xchem-CAR] selected and choose **"Ok"**/**"Open"**
- Your container should start to build, click on the popup notification at the bottom right of visual studio to view the log/progress

### <a name="TimetoLaunch"></a>Time to Launch 

- Open a new terminal that you can interact with. if the terminal is visible at the bottom of the screen click on the plus "create new integrated terminal" or use the keyword shortcut "**Ctrl+Shift+`**" button or use the adjacent "split terminal" (or "**Ctrl+Shift+5**") button to see the new terminal adjacent to the current terminal
- you should now be in the container running Debian Linux
- in the new terminal type (Terminal 1):
  > `cd /container/.devcontainer/` <br> 
  > `chown -R root launch-backend.sh launch-frontend.sh` <br>
  > `./launch-backend.sh` <br>
  > `./launch-frontend.sh` <br>
  
- to start the Vite server for the frontend, in a separate terminal inside your development container:
  >`cd /container/src/frontend/`<br>
  >`npm run dev`<br>z

- launch the django server in a separate terminal:
  >`python3 manage.py runserver`<br>`


for future launches, you will not need to perform the migrations and install the node packages; you will only need to launch the Django server by running: 
  >`cd /container/src/`<br> 
  >`python3 manage.py runserver`<br>

to upload files in CAR, you need to start a Celery worker in a separate terminal inside your development container:

- open a new terminal the same way as last time ([see Time to Launch](#TimeToLaunch))
- in the new terminal type:
  > `cd /container/src/`<br> 
  >`celery -A CAR worker -l info`<br>

if you make any changes to the Django models, you will need to run the the migrations again in the CAR directory:

- makemigrations and migrate the Django models:
  >`python3 manage.py makemigrations` <br> 
  >`python3 manage.py migrate` <br>
  
  or

- use the launch-backend.sh script

# Opening the application

at the end of the step "[Time to Launch](#TimeToLaunch)" an address to use the visual interface should have been displayed ("http://127.0.0.1:3000/"), Ctrl+Click on the link in terminal or copy and paste the link into your web browser to use the CAR interface

> http://127.0.0.1:3000/

the Django server and Rest API can be found at: 

> http://127.0.0.1:8000/

### <a name="Troubleshooting"></a>Troubleshooting

WLS2/Windows users have reported problems with instructions in "Time to Launch":

### <a name="File Permissions"></a>File Permissions

Sometimes, when a new branch is created on Github, try to `git pull`, get a keyerror for accessing the environment variables in the development container or trying to run the `npm run dev` yielding a `Error: EACCES: permission denied`, you need to change the file permissions on the repo folder.

In your WSL2 Linux terminal **NB outside** your dev container:

- change the file permissions of the repository folder
  `sudo chown -R <username> <path to xchem-CAR>`
