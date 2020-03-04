# ppknn
Privacy Preserving k-ish Nearest Neighbors
===

Intruduction
---

This code implements the k-ish Nearest Neighbors (kish NN) classifier as described in PoPETS 2020 (to appear).

Requirements
---
The library requires a Linux environment, a c++11 compiler. In addition it requires the HElib, NTL, GMP and LiPHE libraries.

Installation
---
To compile the code, edit the Makefile and update the directories of of the required libraries.

then run
```
make
```

Databases
---
The code comes with te databases we tested on:
- Breast Cancer Database (also projected onto 3,5 dimensions)
- Car evaluation database


Running
---
To compare the kish NN classifier with the kNN classifier:

```
./test_kish_vs_k --in=databases/cars.csv
```

To test how well k-ish performs with the Gaussian assumption:
```
./test_gaussian_kish --in=databases/cars.csv --res=100
```


To benchmark the FHE k-ish NN classifier on Gaussian data:

```
./test_helib --L=54 --in=databases/breast_cancer_classification_3d.csv --t=16 --res=50
```



The tests were successful on Ubuntu 16 with
- g++ 5.4.0
- HElib
- NTL 10.5.0
- GMP 6.1.0



Feedback / Citation
---
Please send any feedback to Hayim Shaul (<hayim@mit.edu>).
If you use this library please cite it.
If an algorithm you used was described in a paper, please cite that paper as well.

License
---
The software is released under the MIT License as detailed in `LICENSE`.



