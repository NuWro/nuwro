It is assumed that you have already installed ROOT v30 or higher 
If not please visit [cern root](https://root.cern/install/build_from_source/) for building root from source 

To setup ROOTEGPythia6 run `./setupROOTEGPythia6.sh`. In order to do so first make `setupROOTEGPythia6.sh` an executable file. 

```bash
chmod +x setupROOTEGPythia6.sh
```

The chmod +x command changes the permissions of a file to make it executable. The chmod command stands for change mode, and the +x option adds the execute permission for all users.

Then simply run `./setupROOTEGPythia6.sh`

It clones the repo in `${HOME}` directory, and build it in  `ROOTEGPythia6_DIR/build`
It then setup the environmental variable `ROOTEGPythia6_ROOT` in your `.zshrc` (or `.bashrc`) file. 

Once the script ran successfully. Now proceed with installation of nuwro 

The commands to install nuwro is same as before

```bash 
make install 
```
