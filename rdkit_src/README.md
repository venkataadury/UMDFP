## Interfacing with RDKit
Some extra features are only available when interfacing with RDKit. [RDKit is available for C++](https://www.rdkit.org/docs/GettingStartedInC++.html) as well.<br/>
To be able to compile these modules, you will need to install development files for rdkit. This can be done using conda/mamba
```
conda activate my-env
conda install -c conda-forge librdkit-dev eigen libboost-devel
```
