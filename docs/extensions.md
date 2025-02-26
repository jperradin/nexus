Extensions
==========

## Using Extensions

Nexus-CAT supports various extensions to analyze different systems and cluster connectivities. Each extension provides specific settings and methods tailored to the system it is designed for. Here are some of the available extensions:

### SiOz Extension

The `SiOz` extension is used to analyze SiO2 systems. It includes methods to calculate coordination numbers, cluster connectivities, and more.

### SiSi Extension

The `SiSi` extension is designed for analyzing silicon clusters. It provides functionalities to calculate coordination numbers and identify different types of silicon clusters.

### OO Extension

The `OO` extension focuses on oxygen clusters. It includes methods to calculate coordination numbers and identify different types of oxygen clusters.

### Na Extension

The `Na` extension is used for analyzing sodium clusters. It provides methods to calculate coordination numbers and identify sodium-oxygen clusters.

## Writing Your Own Extensions

To write your own extension for Nexus-CAT, follow these steps:

### 1. Create a New Python File

Create a new Python file in the `extensions` directory. Name the file according to the system you are analyzing, for example, `MyExtension.py`.

### 2. Define Supported Elements

List the elements supported by your extension:

```python
LIST_OF_SUPPORTED_ELEMENTS = ["Element1", "Element2"]
```

### 3. Create Atom Subclasses

Define subclasses for each element, inheriting from the `Atom` class:

```python
class Element1(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "Element2"])
```

### 4. Implement Required Methods

Implement the required methods for your extension:

- `transform_into_subclass(atom: Atom) -> object`
- `get_connectivity(cluster_settings) -> list`
- `get_default_settings(criteria="bond") -> dict`

### 5. Register the Extension

Register your extension in the `__init__.py` file:

```python
from . import MyExtension
```

### Example

Here is an example of a simple extension:

```python
LIST_OF_SUPPORTED_ELEMENTS = ["X", "Y"]

class X(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "Y"])

class Y(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "X"])

def transform_into_subclass(atom: Atom) -> object:
    if atom.get_element() == "X":
        return X(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    elif atom.get_element() == "Y":
        return Y(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    else:
        raise ValueError(f"ERROR: Atom {atom.element} - {atom.id} cannot be transformed into X or Y object.")

def get_connectivity(cluster_settings) -> list:
    return ["X-Y"]

def get_default_settings(criteria="bond") -> dict:
    from ..settings.parameter import Parameter, ClusterParameter

    list_of_elements = [
        {"element": "X", "alias": 1, "number": 0},
        {"element": "Y", "alias": 2, "number": 0},
    ]

    dict_cluster_settings = {
        "criteria": criteria,
        "connectivity": ["X", "Y"],
        "polyhedra": [[4, 4], [5, 5], [6, 6]],
    }

    list_of_cutoffs = [
        {"element1": "X", "element2": "Y", "value": 3.00},
    ]

    dict_settings = {
        "extension": Parameter("extension", "MyExtension"),
        "structure": Parameter("structure", list_of_elements),
        "cluster_settings": ClusterParameter("cluster_settings", dict_cluster_settings),
        "cutoffs": Parameter("cutoffs", list_of_cutoffs),
    }

    return dict_settings
```

By following these steps, you can create custom extensions to analyze different systems with Nexus-CAT.
### Tips for Writing Extensions

Creating your own extension for Nexus-CAT can be a straightforward and rewarding process. Here are some tips to help you along the way:

1. **Start with an Existing Extension**: If you're unsure where to begin, look at the code for existing extensions like `SiSi` or `OO` for the most basic one and `SiOz` for a more in-depth analysis. You can copy and modify pieces of code to fit your needs.

2. **Keep It Simple**: Focus on implementing the essential methods first. You can always add more features later.

3. **Test Frequently**: Regularly test your extension with different datasets to ensure it works as expected.

4. **Ask for Help**: If you run into issues, don't hesitate to ask for help from Github issues or refer to the documentation.

By following these tips and the steps outlined above, you'll be able to create powerful extensions to analyze a variety of systems.