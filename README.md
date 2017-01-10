To get the most up to date version of NuWro use:

```
git clone https://github.com/NuWro/nuwro.git
```

To get one of the stable releases use:

```
git clone -b [release tag] https://github.com/NuWro/nuwro.git
```

e.g.

```
git clone -b nuwro_17.01 https://github.com/NuWro/nuwro.git
```

---

1. Install [ROOT](https://root.cern.ch/) with [PYTHIA6](https://pythia6.hepforge.org/)
2. Set ROOTSYS to the folder where ROOT is installed
3. In this directory type:

  ```
  make
  ```

4. Add nuwro/bin to your PATH
5. Copy data/params.txt to your working directory and modify as needed.
6. Type:

  ```
  nuwro
  ```

7. To analyse the output type:

  ```
  myroot eventsout.root
  ```

For more detailed information see [tutorial](http://www.ift.uni.wroc.pl/~tgolan/talks/NuWro_howto.pdf) or [wiki](https://github.com/NuWro/nuwro/wiki)
