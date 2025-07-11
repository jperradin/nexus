from setuptools import setup, find_packages

setup(
    name="nexus-cat",
    version="2.0.0",
    description="Nexus is a Cluster Analysing Toolkit package.",
    author="The Disordered Organization",
    author_email="julien.perradin@protonmail.fr",
    url="https://github.com/TheDisorderedOrganization/nexus",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "tqdm>=4.50.0",
        "colorama>=0.4.4",
        "numba>=0.53.0",
        "psutil>=5.8.0",
    ],
    python_requires=">=3.9",
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
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Operating System :: OS Independent",
    ],
    project_description=open("README.md").read(),
)
