# NuWro

NuWro is a Monte Carlo neutrino event generator, created at the University of Wroclaw. It simulates neutrino-nucleon and neutrino-nucleus reactions for energies from threshold to TeV. The generator has a detector geometry module and can handle realistic neutrino beams, which make it suitable to use in neutrino experiments. [ROOT](https://root.cern.ch/) framework is used to store generated events.

More information can be found in [NuWro User Guide](https://nuwro.github.io/user-guide/).

# Installation

* Install [ROOT](https://root.cern.ch/) with [PYTHIA6](https://pythia6.hepforge.org/)
* Set *ROOTSYS* to the folder where ROOT is installed
* In NuWro directory type:

  ```
  make
  ```

For more details see [user-guide/getting-started/installation](https://nuwro.github.io/user-guide/getting-started/installation/).

# Usage

* Add *nuwro/bin* to your *PATH*
* Copy *data/params.txt* to your working directory and modify as needed (see [user-guide/getting-started/parameters](https://nuwro.github.io/user-guide/getting-started/parameters/))
* Type

  ```
  nuwro
  ```

* To analyse the output type

  ```
  myroot eventsout.root
  ```

For more details see [user-guide/getting-started/running](https://nuwro.github.io/user-guide/getting-started/running/) and [user-guide/analyze-output](https://nuwro.github.io/user-guide/analyze-output/basics/).

# Releases

Every stable release is tagged using the year followed by the month of the release, e.g. nuwro_17.01 is a stable release from January 2017. Not-tagged versions should be treated as release candidates.

# Credits

Currently involved in the project:

* Jan Sobczyk
* Cezary Juszczak
* [Tomasz Golan](http://www.ift.uni.wroc.pl/~tgolan/)
* [Kajetan Niewczas](http://www.ift.uni.wroc.pl/~kniewczas/)
* [Krzysztof Graczyk](http://www.ift.uni.wroc.pl/~kgraczyk/)

Former members with significant contribution:

* [Jarosław Nowak](http://www.lancaster.ac.uk/physics/about-us/people/jaroslaw-nowak)
* Jakub Żmuda

and

* Artur Kobyliński
* Maciej Tabiszewski
