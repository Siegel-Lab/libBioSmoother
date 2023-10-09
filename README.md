# libSmoother - The C++ library behind Smoother


[Smoother](https://www.github.com/Siegel-Lab/smoother "Smoother GitHub") is an interactive analysis and visualization software for contact mapping data.
This repository, libSmoother, contains the backend C++ library that powers Smoother.
The library takes care of all the time-sensitive computations.
It is written in C++17 and offers seamless integration with Python 3.

## Quick Start

create & activate a new environment (optional)
```
conda create -y -n smoother python=3.9
conda activate smoother
```

Install libSmoother (and all requirements) using pip. LibSmoother runs under Windows, Linux, and MacOS.
```
pip install libbiosmoother
```

Download 2 example smoother indices.
```
wget https://syncandshare.lrz.de/getlink/fi4kLPLjRjMTbRnij7PtyB/t_brucei_hi_c.smoother_index.zip
#wget https://syncandshare.lrz.de/getlink/fiMo5Zsj8baXjXpzD8Whic/m_musculus_radicl_seq.smoother_index.zip

conda install unzip
unzip t_brucei_hi_c.smoother_index.zip
#unzip m_musculus_radicl_seq.smoother_index.zip
```

Export a picture from one of the indices
```
libbiosmoother export t_brucei_hi_c-c -f png
#libbiosmoother export m_musculus_radicl_seq -f png
```

## Full Documentation

Since [Smoother](https://www.github.com/Siegel-Lab/smoother "Smoother GitHub") and libSmoother are designed as one package, they share documentation.
Hence, for more information and in-depth instructions, check out the [Smoother's manual](https://biosmoother.readthedocs.io/ "Smoother's Manual").

## Cite

If you use Smoother in your research, please cite:
@todo