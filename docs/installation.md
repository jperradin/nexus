Installation
============

To install Nexus-CAT, follow these steps:

### Prerequisites

- Python 3.9 or higher
- Required dependencies: `numpy`, `scipy`, `tqdm`, `natsort`

### Installation Steps

#### Basic Installation

To install Nexus-CAT as a package, you can use pip:

```bash
pip install nexus-cat
```

Note: the package does not auto upgrade itself, please run the following command to upgrade to the latest version:

```bash
pip install nexus-cat --upgrade
```

#### Installation from the source code

If you want to install the package from the source code to implement your extensions for example, you can clone the repository:

```bash
git clone https://github.com:TheDisorderedOrganization/nexus.git
```

Then install the package in development mode:

```bash
cd nexus
pip install -e . # will use the setup.py file to install the package
```

You have successfully installed Nexus-CAT on your system. For more details, refer to the [Getting Started](getting_started.rst) guide.