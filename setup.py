from setuptools import setup, find_packages

setup(
    name="nexus-cat",
    version="1.0.0",
    description="Nexus is a Cluster Analysing Toolkit package for atomic systems.",
    author="Julien Perradin",
    author_email="julien.perradin@protonmail.fr",
    url="https://github.com/jperradin/nexus",
    packages=find_packages(),
    install_requires=["numpy", "scipy", "tqdm", "natsort"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Materials Science",
        "Topic :: Scientific/Engineering :: Computational Physics",
        "Topic :: Scientific/Engineering :: Computational Chemistry",
        "Topic :: Scientific/Engineering :: Computational Materials Science",
        "Topic :: Scientific/Engineering :: Computational Physics",
        "Topic :: Scientific/Engineering :: Computational Chemistry",
        
    ],
    project_description=open("README.md").read(),
)
