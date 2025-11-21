🛑 It is assumed that you have already installed `ROOT v30` or higher  
If not please visit [cern root](https://root.cern/install/build_from_source/) for building root from source 

💡 You have to run this script only once setup ROOTEGPythia6.

## How to execute the script

First make `setupROOTEGPythia6.sh` an executable file. 
```bash
chmod +x setupROOTEGPythia6.sh
```

Then

run 
```bash
./setupROOTEGPythia6.sh
```
in your terminal. The chmod +x command changes the permissions of a file to make it executable. The chmod command stands for change mode, and the +x option adds the execute permission for all users.

## Details of the script

The script clones this repository [`ROOTEGPythia6`](https://github.com/luketpickering/ROOTEGPythia6.git) 

🎁 many thanks to [Luke Pickering](https://www.dunescience.org/facesofdune/luke-pickering/)

For more details please visit [`ROOTEGPythia6`](https://github.com/luketpickering/ROOTEGPythia6.git).

By default, it clones and build the repo in `${HOME}` directory, however one can also provide an absolute path where the repo should be cloned and build

💡 **For convenience, it is suggested to clone and build this repo outside `nuwro/` directory as deleting the `nuwro/` won't affect the `ROOTEGPYthia6`**

The build material is in `/path/to/ROOTEGPythia6/build` directory

It then setup the environmental variable `ROOTEGPythia6_ROOT` in your `.zshrc` (or `.bashrc`) file (depending on your shell ...)

Once the script ran successfully. Now proceed with installation of nuwro 


📝: Please note that the commands to install nuwro is same as before

```bash 
make install 
```

If you face any issues with this script, please contact the following :

**Hemant Prasad**, [Email](hemant.prasad.pl@gmail.com) 📬, Faculty of Physics and Astronomy, Plac Maxa Borna 9, University of Wrocław, Poland

