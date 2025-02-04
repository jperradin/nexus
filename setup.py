from setuptools import setup, find_packages

setup(
    name='nexus-cat',
    version='0.1.11',
    description='Nexus is a Cluster Analysing Toolkit package for atomic systems.',
    author='Julien Perradin',
    author_email='julien.perradin@umontpellier.fr',
    url='https://github.com/JulienPerradin/nexus',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'tqdm',
        'natsort'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    project_description=open('README.md').read(),
)
