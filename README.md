[![CI](https://github.com/Waztom/xchem-CAR/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/Waztom/xchem-CAR/actions/workflows/ci.yml)

[![Docker Hub](https://img.shields.io/docker/v/waztom/car?label=Docker%20Hub&logo=docker&sort=date)](https://hub.docker.com/r/waztom/car)
[![Docker Pulls](https://img.shields.io/docker/pulls/waztom/car?logo=docker)](https://hub.docker.com/r/waztom/car)
[![Docker Image Size](https://img.shields.io/docker/image-size/waztom/car?sort=date&logo=docker)](https://hub.docker.com/r/waztom/car)

[![Python 3.15](https://img.shields.io/badge/python-3.15-blue?logo=python&logoColor=white)](https://python.org)
[![Django 6.1](https://img.shields.io/badge/django-6.1-green?logo=django&logoColor=white)](https://djangoproject.com)
[![React 18](https://img.shields.io/badge/react-18-61DAFB?logo=react&logoColor=white)](https://react.dev)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow)](LICENSE)

The containers use Python 3.15 with pip-managed dependencies in an isolated
`/opt/venv`, including RDKit and psycopg2. Rebuild existing devcontainers to
apply the change. See the [Python environment notes](docs/source/reference/contributing.rst).

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

**Upgrading an existing installation from PostgreSQL 10?** Choose one of the
[PostgreSQL 17 options below](#postgresql-10-to-17) before launching CAR.

- start Visual Studio remote container with **Ctrl + Shift + P** and type **"Remote-containers: Open folder in container"** then click on that option. <br> ensure you have the repository folder [your file path/xchem-CAR] selected and choose **"Ok"**/**"Open"**
- Your container should start to build, click on the popup notification at the bottom right of visual studio to view the log/progress

### PostgreSQL 10 to 17

The development container now uses **PostgreSQL 17 (PG17)** instead of
**PostgreSQL 10 (PG10)**. A Docker *volume* stores database data separately
from the container, so rebuilding a container does not upgrade its database.
PG17 cannot read a PG10 volume directly.

Choose **one** option:

| Option | Choose this if... | Result |
| --- | --- | --- |
| [1. Start fresh](#option-1-start-fresh-with-pg17) | You do not need your old local database records. | CAR starts with an empty PG17 database. |
| [2. Keep your data](#option-2-migrate-your-pg10-database-to-pg17) | You need your existing users, experiments and other database records. | Export the PG10 database and restore it into PG17. |

The updated [development Compose configuration](.devcontainer/docker-compose-dev.yml)
uses a new volume named `pgdata17` instead of `pgdata`. Docker usually prefixes
these names with a project name, for example `xchem-car_devcontainer_pgdata`.
The old volume is **not deleted or automatically migrated**.

**Safety first**

- These instructions are for the **local development container**, not a shared
  or production database.
- Do not delete the old volume, run `docker compose down -v`, or use Docker
  Desktop's volume cleanup/prune options. Keep the old data until you have
  checked the new installation and saved a backup.
- Never point PG17 at the old `pgdata` volume or simply rename it to `pgdata17`.
- A database backup does **not** include uploaded files in `src/media` or other
  files stored outside PostgreSQL. Back these up separately and keep the same
  checkout when migrating.

#### Option 1: Start fresh with PG17

This leaves your PG10 volume untouched, but **its records will not appear in
the new database**. Old database users will not be available either.

1. Stop your running Django server and Celery workers. If VS Code is inside
   the container, use **Ctrl+Shift+P -> Dev Containers: Reopen Folder Locally**.
2. Make sure your checkout includes the updated
   [development Compose configuration](.devcontainer/docker-compose-dev.yml):
   the database image should be `postgres:17-alpine` and its volume `pgdata17`.
   Unlock the secrets as described [above](#UnlockingSecrets).
3. Use **Ctrl+Shift+P -> Dev Containers: Rebuild and Reopen in Container**
   (or **Dev Containers: Rebuild Container** if already connected). Wait for
   the build to finish.
4. Continue with [Time to Launch](#TimetoLaunch). The backend launch script
   runs Django migrations to create the tables in the new database.
5. If you need a Django administrator account, open a terminal **inside** the
   devcontainer and run:

   ```bash
   cd /container/src
   python3 manage.py createsuperuser
   ```

If you already have a `pgdata17` volume, rebuilding **reuses its data**; it
does not empty it. If that volume contains data you want to replace, stop here
and ask the maintainer to help identify and back up the correct volume before
resetting it. Do not delete volumes by guesswork.

#### Option 2: Migrate your PG10 database to PG17

Allow time for a backup and restore; larger databases take longer. Do not use
CAR or allow background workers to write data during this process.
**Stop at any error** and keep the error message for the maintainer; do not
continue with a failed or incomplete backup/restore.

Run the Docker commands below in a terminal **on your computer, outside the
devcontainer**, from the repository root (the folder containing this README).
On Windows, use a **WSL terminal**, not PowerShell or Command Prompt. Docker
Desktop/Engine must be running, and the repository secrets must be unlocked.
Replace names beginning with `YOUR_` with the actual names from your machine.

**1. Identify and stop the old database**

- Stop Django and all Celery workers, then reopen the folder locally in VS Code
  using **Dev Containers: Reopen Folder Locally**.
- In Docker Desktop, find this repository's old database container. Its image
  should be PostgreSQL 10. The following command also lists database containers
  from Compose projects:

  ```bash
  docker ps -a --filter label=com.docker.compose.service=db --format 'table {{.Names}}\t{{.Image}}'
  ```

- Find the old container's database volume in Docker Desktop's **Inspect**
  view, or run:

  ```bash
  docker inspect YOUR_OLD_DB_CONTAINER --format '{{range .Mounts}}{{println .Name .Destination}}{{end}}'
  ```

  Record the full volume name mounted at `/var/lib/postgresql/data` (a trailing
  `/` is also normal). It normally ends in `_pgdata`, **not** `_pgdata17`.
  If the old container has already been removed, use Docker Desktop's
  **Volumes** view or `docker volume ls` to find the old project's volume.
  If you are unsure which volume belongs to CAR, ask the maintainer before
  proceeding.
- Stop the old database container, if it still exists:

  ```bash
  docker stop YOUR_OLD_DB_CONTAINER
  ```

  Also stop this project's other application/worker containers in Docker
  Desktop. Do not stop unrelated projects, and never run two database
  containers against the same volume.

**2. Open the old volume with PG10 and make a backup**

First check that the volume exists. If this command reports "no such volume",
stop and correct the name:

```bash
docker volume inspect YOUR_OLD_PG10_VOLUME
```

Start a temporary PG10 container using that exact volume name:

```bash
docker run -d --name car-pg10-migration --network none \
  --env-file .devcontainer/.devcontainerenv \
  --mount type=volume,source=YOUR_OLD_PG10_VOLUME,target=/var/lib/postgresql/data \
  postgres:10-alpine
docker logs car-pg10-migration
```

This publishes no network ports. Wait until the logs say the database is
ready to accept connections; repeat the logs command if necessary. If they
report an incompatible database version, stop the temporary container and ask
for help. Do not initialize or erase the volume. If the temporary container
name is already in use, inspect that container before proceeding.

The commands use the database name and user from the unlocked environment
file (or PostgreSQL's defaults). These must match the **old** database;
changing the environment file does not change existing database credentials.
If you previously used custom database names, roles or extensions, ask the
maintainer to adapt these steps.

Create a new, private backup folder outside the repository, then export CAR's
database using the **PG17 client tool**. This client can read PG10 and prepares
the dump for the newer server. It connects only to the temporary container's
private network namespace; no database port is exposed on your computer.

```bash
mkdir -m 700 "$HOME/car-pg10-migration"
```

If the folder already exists, stop and choose a new folder name for all
subsequent commands; do not overwrite an earlier backup. After the folder
has been created successfully, run the export:

```bash
docker run --rm --network container:car-pg10-migration \
  --env-file .devcontainer/.devcontainerenv \
  postgres:17-alpine sh -c '
  db_user="${POSTGRES_USER:-postgres}"
  export PGPASSWORD="$POSTGRES_PASSWORD"
  exec pg_dump -h 127.0.0.1 -U "$db_user" -d "${POSTGRES_DB:-$db_user}" --format=custom
' > "$HOME/car-pg10-migration/car.dump"
```

Check the export succeeds before continuing. Keep the dump private: it
contains application data and password hashes.

Check that PG17's restore tool can read the backup's table of contents:

```bash
docker run --rm --network none \
  --mount "type=bind,source=$HOME/car-pg10-migration,target=/backup,readonly" \
  postgres:17-alpine pg_restore --list /backup/car.dump
```

Expect a list of database objects, not an error. This checks the archive
format; the restore and application checks below are still essential.
Once the export and archive check succeed, stop the temporary container:

```bash
docker stop car-pg10-migration
```

**3. Build PG17, but do not launch CAR yet**

Follow steps 2 and 3 of [Option 1](#option-1-start-fresh-with-pg17) to rebuild
with the updated configuration. **Do not run the backend launch script,
Django migrations, the Django server or Celery yet.** The restore needs a
new database without CAR tables.

In your **host/WSL terminal**, list the database containers again:

```bash
docker ps --filter label=com.docker.compose.service=db --format 'table {{.Names}}\t{{.Image}}'
```

Identify this checkout's `postgres:17-alpine` container. Confirm in Docker
Desktop or with the `docker inspect` command from step 1 that it uses the
new volume ending in `_pgdata17`, not the old volume. Check readiness:

```bash
docker logs YOUR_PG17_DB_CONTAINER
docker exec YOUR_PG17_DB_CONTAINER sh -c '
  db_user="${POSTGRES_USER:-postgres}"
  exec psql -U "$db_user" -d "${POSTGRES_DB:-$db_user}" -c "\dt"
'
```

Wait for the logs to say the database is ready. The table check should report
that no relations were found. If CAR tables already exist, **stop**: do not
restore over them or add `--clean`. Ask the maintainer to help prepare a
separate empty target database without losing existing data.

**4. Restore the backup into PG17**

Still in your **host/WSL terminal**, run:

```bash
docker exec -i YOUR_PG17_DB_CONTAINER sh -c '
  db_user="${POSTGRES_USER:-postgres}"
  exec pg_restore -U "$db_user" -d "${POSTGRES_DB:-$db_user}" \
    --no-owner --no-acl --exit-on-error
' < "$HOME/car-pg10-migration/car.dump"
```

A successful restore normally prints nothing. Wait for the command to finish.
This restores one application database, assigning objects to the configured
database user; it does not migrate cluster-wide roles or old permission grants.
If it fails, do not launch CAR or retry into the partly restored database:
save the error and ask the maintainer for help with a fresh target.

**5. Apply application updates and check your data**

- Continue with [Time to Launch](#TimetoLaunch) **inside the devcontainer**.
  The backend launch script applies any newer Django migrations to the restored
  database. This is separate from the PostgreSQL version upgrade.
- Sign in with an existing CAR account and check familiar experiments,
  records and uploaded files. Confirm the expected data is present and that
  normal application operations work before treating the migration as complete.
- Keep the old PG10 volume and your dump until these checks pass and you have
  a separate safe copy of the backup. You can then remove **only the temporary
  container** (this command does not delete the named database volume):

  ```bash
  docker rm car-pg10-migration
  ```

If you need to return to PG10, stop the new application and workers and ask
the maintainer to help run the previous application version against the
preserved PG10 volume. Do not attach PG10 to `pgdata17`. Changes made in PG17
after migration are not copied back to PG10.

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
