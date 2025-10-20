The files within this folder are hosted on Google drive

Running the tests requires the google drive folder be placed here, otherwise data will 
be missing

The folder can be found at:

https://drive.google.com/drive/folders/1Zm0fjwL2bZPdUGRAdEJeaNmtiklNI8RV?usp=sharing


There are two options to get this data:

1. Download the data, unzip the folder, and place the contents here.

2. Use `rclone` to copy the data. This is a little more involved, but can be faster in
the long-run if you are going to be doing development and/or updating the data within
this folder. First, make sure the folder is appearing in "Shared with me" on Google Drive.
  - Run `rclone config` if you haven't already. Follow the interactive process.
  - Either run `rclone config file` and add the line `shared_with_me = true` to the 
    section for Google Drive, or add the `--drive-shared-with-me` flag to each of the
    listed commands.
  - The rclone command to download is `rclone copy [REMOTE]:PyITMData/ .`. Perform a
    dry-run or run interactively if you are unsure.
  - To upload run:
    `rclone copy --exclude='*.csv' --exclude='*.md' . [REMOTE]:PyITMData/` (make sure
    to change to your remote's name & double-check the paths). You may with to run 
    interactively or do a dry run first.


--------------------

After downloading the data, run the script "./config.sh" from this directory.

```shell
chmod +x config.sh
./config.sh
```

The script will only configure the paths in the lookup file. It runs sed (Stream EDitor)
to set the paths in the template satmap file to the absolute path of the files.

