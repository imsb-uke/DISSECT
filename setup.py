from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='dissect_deconv',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.19.5',
        'scikit-learn>=1.1.2',
        'pandas>=1.4.3',
        'anndata>=0.8.0',
        'h5py>=3.6.0',
        'scanpy>=1.9.1',
        'shap>=0.40.0',
        'Cython>=0.29.28',
        'matplotlib>=3.5.0',
        'setuptools>=59.4.0',
        'tqdm>=4.62.3',
#        'keras==2.7.0',
#        'tensorflow>=2.7.0',
    ],
    ext_modules=cythonize("dissect_deconv/PropsSimulator/*pyx"),
    keywords=[
        "bioinformatics",
        "deep learning",
        "machine learning",
        "single cell sequencing",
        "bulk rna sequencing",
        "spatial transcriptomics",
        "deconvolution",
    ],
    python_requires=">=3.8.6",
    author="Robin Khatri",
    author_email="robin.khatri@zmnh.uni-hamburg.de",
    url="https://github.com/imsb-uke/DISSECT",
    
    license="MIT License",
)
